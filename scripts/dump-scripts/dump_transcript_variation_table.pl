#!/usr/local/bin/perl 

=head1 NAME

dump_transcript_variation_table.pl - generate variation table like the one on the browser gene page 
	

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
use Getopt::Long;
use Pod::Usage;


=head1 SYNOPSIS

dump_transcripts.pl  [options] 
 
 Options:
    --help		help message
    --man		full documentation
    --registry          the registry file for database connections
    --speceis 		which species to dump
    --id_file           the file containing transcript ids
    --debug		debug
    --format		two types, trp or var
 
=head1 OPTIONS

=over 4


=item B<--exclude>

    Genes to ignore.  


=item B<--exclude-clone>
    
    Clones to ignore 
    This may be a comma-separated list.
    
=item B<--registry>

The registry file for ensembl databases

=item B<--species> 

supply the species name whose transcripts are to be dumped


=item B<--id_file>

the file containing transcript ids

=item B<--format>

reporting format
var:VariantID	Chr:bp	vf_allele	Alleles	ConseqType	Transcript
trp:Transcript(strand)	Allele(transcript_allele)	ConseqType	PositionInTranscript	PositionInCDS	PositionInProt	AminoAcid	Codons
=item B<--help> 

print a help message and exit

=item B<--man> 

print documentation and exit

    
=back


=cut


# Osmh63.01G000010_02.1 as example
# format 1
#Transcript (strand)	Allele (transcript allele)	Consequence Type	Position in transcript	Position in CDS	Position in protein	Amino acid	Codons
#format 2
#Variant ID	Chr: bp	vf_allele	Alleles	Conseq. Type	Transcript

my ($species, $registry, $idfile, $debug, $format);

{  #Argument Processing
  my $help=0;
  my $man=0;

  GetOptions( "help|?"=>\$help,"man"=>\$man
	      ,"species=s"=>\$species
	      ,"registry=s"=>\$registry
	      ,"id_file=s"=>\$idfile
	      ,"debug"=>\$debug
	      ,"format=s"=>\$format
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

our $trpt_adaptor = $CORE_DBA->get_TranscriptAdaptor;
our $slice_adaptor = $CORE_DBA->get_SliceAdaptor; 
our $tv_adaptor = $VAR_DBA->get_TranscriptVariationAdaptor;
our $vf_adaptor = $VAR_DBA->get_VariationFeatureAdaptor;

my @print_out;

if( $format eq 'trp'){
	
	my @transcripts;
	if( $idfile ){
		open my $fh, $idfile or die "cannot open file $idfile";
		@transcripts = map{ $trpt_adaptor->fetch_by_stable_id( $_ ) }
	 		grep { /\S+/} map { chomp; s/\s+//g; $_;} <$fh>;
		close $fh;
	}else{
		my $slices = $slice_adaptor->fetch_all('toplevel');
		map{ push @transcripts, @{$_->get_all_Transcripts} } @$slices;
	}

	my %count;

	foreach my $t (@transcripts) {
  		my $id = $t->stable_id;
  		print "! $id\n" if $debug;
  		push @print_out, map{join "\t", @_} @{get_transcript_variations($t)};  
  		$count{trpt}++;
  		last if $debug;
	}

	for my $k (sort keys %count){
  		print "$k = $count{$k}\n";
	}

}elsif( $format eq 'var' ){

	my $slices = $slice_adaptor->fetch_all('toplevel');
        #map{ push @transcripts, @{$_->get_all_Transcripts} } @$slices;
	
	for my $slice(@$slices){
		push @print_out, map{join "\t", @_} @{get_variation_features($slice)};
		last if $debug;
	}
}else{

	warn("Need to choose format, var or trp\n");
}

for ( @print_out){
	print "$_\n";
}

sub get_transcript_variations {
  	my $transcript   = shift;

	my @tvs = @{ $tv_adaptor->fetch_all_by_Transcripts([$transcript]) };
	my $tvs_cnt = scalar @tvs;
	print "Get $tvs_cnt transcript variation pairs for " . $transcript->stable_id ."\n" if $debug;

	my @tvs_report;
	for my $tv ( @tvs ){
                my $vf = $tv->variation_feature();
        
                push @tvs_report,[
				$tv->transcript_stable_id,
                                $vf->variation_name,
                                $vf->seq_region_name . '(' .$vf->seq_region_strand .'):'.
                                	$vf->seq_region_start . '-' . $vf->seq_region_end,  
				$vf->allele_string,
				$vf->alt_alleles,
                                $tv->consequence_type,  
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
                                $_->consequence_type,
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

