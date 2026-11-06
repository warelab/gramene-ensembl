#!/usr/local/bin/perl 

=head1 NAME

restore_sift_scores.pl - populate sift scores in transcript_variation table from protein_function_predictions tables. 
	

=cut

BEGIN {
    $ENV{'GrameneDir'} = '/usr/local/ensembl-live/gramene-live/'; 
    $ENV{'GrameneEnsemblDir'} = '/usr/local/ensembl-live/'; 
}

use lib map { $ENV{'GrameneDir'}."/$_" } qw ( lib/perl );
use lib map { $ENV{'GrameneEnsemblDir'}."/$_" } 
        qw ( bioperl-live
	     ensembl/modules 
	     ensembl-variation/modules
	     ensembl-compara/modules );

use strict;
use warnings;


use Bio::EnsEMBL::DBLoader;
use Bio::EnsEMBL::Registry;
use List::Util qw(max min);
use Bio::EnsEMBL::Variation::Variation;
use Bio::EnsEMBL::Variation::VariationFeature;
use Bio::EnsEMBL::Variation::Utils::Config qw(%ATTRIBS);
use Bio::EnsEMBL::Variation::Utils::Constants qw(%VARIATION_CLASSES);
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw($UPSTREAM_DISTANCE $DOWNSTREAM_DISTANCE);
use Bio::EnsEMBL::Variation::Utils::VariationEffect;
use Bio::EnsEMBL::Variation::TranscriptVariationAllele;
use Bio::EnsEMBL::Variation::DBSQL::ProteinFunctionPredictionMatrixAdaptor;
use Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix;
use Getopt::Long;
use Pod::Usage;
use Digest::MD5 qw(md5_hex);

=head1 SYNOPSIS

restore_sift_scores.pl  [options] 
 
 Options:
    --help		help message
    --man		full documentation
    --registry          the registry file for database connections
    --speceis 		which species to dump
    --debug		debug
 
=head1 OPTIONS

=over 4

=item B<--registry>

The registry file for ensembl databases

=item B<--species> 

supply the species name whose transcripts are to be dumped


=item B<--help> 

print a help message and exit

=item B<--man> 

print documentation and exit

    
=back


=cut


my ($species, $registry, $debug);

{  #Argument Processing
  my $help=0;
  my $man=0;

  GetOptions( "help|?"=>\$help,"man"=>\$man
	      ,"species=s"=>\$species
	      ,"registry=s"=>\$registry
	      ,"debug"=>\$debug
	    )
    or pod2usage(2);
  pod2usage(-verbose => 2) if $man;
  pod2usage(1) if $help;
}

# Load the ensembl file
Bio::EnsEMBL::Registry->load_all( $registry );
our $CORE_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'core' );
$CORE_DBA || pod2usage( "\nNo core DB for $species set in $registry\n" );

our $VAR_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'variation' );
$VAR_DBA || pod2usage( "\nNo core DB for $species set in $registry\n" );

our $tadaptor = $CORE_DBA->get_TranscriptAdaptor;
our $slice_adaptor = $CORE_DBA->get_SliceAdaptor; 
our $tvadaptor = $VAR_DBA->get_TranscriptVariationAdaptor;
our $vf_adaptor = $VAR_DBA->get_VariationFeatureAdaptor;
our $v_adaptor = $VAR_DBA->get_VariationAdaptor;
our $pfpm_adaptor = $VAR_DBA->get_ProteinFunctionPredictionMatrixAdaptor;

  # fetch all TranscriptVariations related to a Transcript

#my @test_trpts = ('SORBI_3008G072400.1', 'SORBI_3008G134800.1');
my @test_trpts = ('EES15821', 'EES16183');
#for my $transcript (@{ $tadaptor->fetch_all_by_biotype('protein_coding') }){
			#my $transcript = $tadaptor->fetch_by_stable_id('SORBI_3008G072400.1');
			# SORBI_3008G134800.1

my $transcript;
for my $trpt_sid ( @test_trpts ){
	
	eval{ 
		$transcript = $tadaptor->fetch_by_stable_id($trpt_sid);
	};
 
	if ( $@ ){
		warn("No transcript found for $trpt_sid, $@");
		next;
	}

	unless( defined $transcript){
		warn("No transcript found for $trpt_sid");
                next;
	}	
	my $translation = $transcript->translation;
	my $sequence = $translation->seq;
	#my $translation_id = $translation->stable_id;
			#print "prot for $translation_id is $sequence\n";
	my $md5 = md5_hex(uc($sequence));
					#3af50c6abad1611232170c88581f8943
	
	my ($pfpm, $calc_length) ;
	$pfpm = $pfpm_adaptor->fetch_sift_predictions_by_translation_md5( $md5 );
	unless( defined $pfpm ){
		warn("No pfpm found for $trpt_sid, skip\n$sequence\n$md5\n");
		next;
	}
        $calc_length = &calc_length_pfpm( $pfpm );
	unless( $calc_length ){
		warn ("Failed to get calc_length for trpt_id $trpt_sid\n$sequence\n$md5\n");
		next;
	}
	
	warn( "$trpt_sid, calc_len=$calc_length\n$sequence\n$md5\n");
	for my $tv (@{ $tvadaptor->fetch_all_by_Transcripts( [$transcript]) }){
  		#foreach my $allele (keys %{$tv->hgvs_protein}) {
    		#my $hgvs_notation = $tv->hgvs_protein->{$allele} || 'hgvs notation is NA';
    		#print "$allele $hgvs_notation\n";
  		#}

        	for my $tva (@{ $tv->get_all_alternate_TranscriptVariationAlleles }) {

		        my $missense_variant_flag = 0;
                	my $conseqs = join ",", map { $_->SO_term } @{ $tva->get_all_OverlapConsequences };
                	$missense_variant_flag = 1 if $conseqs =~ /missense_variant/i;
                	next unless $missense_variant_flag;
                                       
                	my $tl_start = $tv->translation_start;
                	#print "translation start is $tl_start\n";
                	warn "$tl_start > $calc_length, skip" and next if ($tl_start > $calc_length);

                	my $aa_change = $tva->pep_allele_string;
                	my ($ref_aa, $alt_aa) = split '/', $aa_change;
                	#print "aa change is $ref_aa -> $alt_aa, at pos $tl_start\n";
                	#print "sequence with respect to the transcript: ", $tva->feature_seq, "\n";
                	#print "sequence with respect to the variation feature: ", $tva->variation_feature_seq, "\n";

                	my ($sift_pred, $sift_score);
                	eval{
                        	($sift_pred, $sift_score) = $pfpm->get_prediction($tl_start, $alt_aa);
                	};
                	if( $@ ){
                        	warn ("Skip: $tl_start, $alt_aa\n $@\n");
                        	next;
                	}

                	print "$trpt_sid protein $tl_start ($calc_length) $ref_aa -> $alt_aa: SIFT pred/score: $sift_pred / $sift_score\n";
      		}

  	} #end of tv

	#last;

} #end of transcript


sub calc_length_pfpm {

	my $pfpm = shift;
	my $debug = shift;

	$pfpm->expand_matrix(); 
        my $binary_string = $pfpm->{matrix};
                
        if (defined $binary_string) {
                my $header_len = 3; # Length of 'VEP'
                my $bytes_per_pos = 40; # 20 amino acids * 2 bytes each

                my $calc_length = (length($binary_string) - $header_len) / $bytes_per_pos;
                warn "Calculated Peptide Length: $calc_length\n" if $debug;
		return $calc_length;
        } else {
                warn "No matrix data found!\n" if $debug;
		return undef;
        }

}


__END__


sub get_transcript_variations {
  	my $transcript   = shift;

	my @tvs = @{ $tv_adaptor->fetch_all_by_Transcripts([$transcript]) };
	my $tvs_cnt = scalar @tvs;
	print "Get $tvs_cnt transcript variation pairs for " . $transcript->stable_id ."\n" if $debug;

	my @tvs_report;
	for my $tv ( @tvs ){
                my $vf = $tv->variation_feature();
        
#Osmh63.01G000010_02     Chr01_11742_C_G 1(1):11742-11742        C/G     ARRAY(0x4343f60)        ARRAY(0x4340df0)        missense_variant 790     415     139     P/A     Cca/Gca       
	        push @tvs_report,[
				$tv->transcript_stable_id,
                                $vf->variation_name,
                                $vf->seq_region_name . '(' .$vf->seq_region_strand .'):'.
                                	$vf->seq_region_start . '-' . $vf->seq_region_end,  
				$vf->allele_string,
				(join "|", @{$vf->alt_alleles}),
                                (join "|", @{$tv->consequence_type}),  
                                $tv->display_consequence,
                                $tv->cdna_start,
				$tv->cds_start,
				$tv->translation_start, 
                                $tv->pep_allele_string, 
                                $tv->codons,
                                ];  
        }
     

  	return \@tvs_report;
}

sub get_variation_features {

     	my $slice = shift;
     	my @variation_features = @{ $vf_adaptor->fetch_all_by_Slice($slice) } ;

	#Variant ID     Chr: bp vf_allele       Alleles Conseq. Type    Transcript
	my @snp_report ;

	for my $vf ( @variation_features ){
		my $tvs = $vf->get_all_TranscriptVariations;

		 my @asnp_tvs =  map{
				[$vf->variation_name,
                                $vf->seq_region_name,
                                $vf->seq_region_start,
                                $vf->seq_region_end,
                                $vf->seq_region_strand,
                                #(join "|", @{$_->consequence_type}),
                                $_->display_consequence,
				$_->transcript_stable_id,
				$_->pep_allele_string, 
				$_->translation_start,
				]
				} @{$tvs};
		push @snp_report, @asnp_tvs;				
	}

     	return \@snp_report;
}


  
########################## subroutines ######################################

__END__


=head1 OUTPUT


=head1 AUTHOR

   Sharon Wei <weix@cshl.edu>

   Gramene Project (www.gramene.org)
   Ware Lab
   Cold Spring Harbor Laboratory

   This is free software and may be distributed on the same terms as Perl itself.

=cut

