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
my $transcript = $tadaptor->fetch_by_stable_id('SORBI_3008G072400.1');
my $sequence = $transcript->translation->seq;
print "prot for SORBI_3008G072400.1 is $sequence\n";
my $md5 = md5_hex(uc($sequence));
print "Current Sequence MD5: $md5\n";
#3af50c6abad1611232170c88581f8943
#my $md5 = '3af50c6abad1611232170c88581f8943';
my $pfpm ;
$pfpm = $pfpm_adaptor->fetch_sift_predictions_by_translation_md5( $md5 );
	#my $tl_start = $tv->translation_start;
	#print "translation start is $tl_start\n";
my ($prediction, $score) = $pfpm->get_prediction(2, 'C');
print "A mutation to 'C' at position 2 is predicted to be $prediction, score $score\n";

# 1. Ensure the matrix is expanded (uncompressed)
# The length of a zipped string won't help with the math.
$pfpm->expand_matrix();

# 2. Get the raw binary string
my $binary_string = $pfpm->{matrix};

# 3. Apply the math
if (defined $binary_string) {
    my $header_len = 3; # Length of 'VEP'
    my $bytes_per_pos = 40; # 20 amino acids * 2 bytes each
    
    my $calc_length = (length($binary_string) - $header_len) / $bytes_per_pos;
    print "Calculated Peptide Length: $calc_length\n";
    
    # 4. Optional: Set it so the getter works next time
    $pfpm->peptide_length($calc_length);
} else {
    print "No matrix data found!\n";
}

print "Matrix Peptide Length: " . $pfpm->peptide_length() . "\n";
print "Matrix Object MD5:    " . $pfpm->translation_md5() . "\n";

print "Querying Position: 472\n";
($prediction, $score) = $pfpm->get_prediction(472, 'E');
print "A mutation to 'E' at position 472 is predicted to be $prediction, score $score\n";

__END__

#my $sequence = "MVICCQLQRLARKTAKIGLIAINSFKKKALHLMGCLCSKGAKDDANATSGRRTPSRKSDS
AADAVSNNGGTAVLNAKAKEKLSGGEKVAVALDARISSGNNAELKGLSGEHVVAGWPSWL
INVAPKAVEGWLPRRADSFEKLAKIGQGTYSIVYKARDLESGKIVALKKVRFVNMDPESV
RFMAREIHILRRLDHPNVIKLEGIVTSRVSQNLYLVFEYMEHDLAGLVATPGLKLTEPQI
KCFVQQLLHGLDHCHKNGVLHRDIKGANLLIDSNGMLKIGDFGLAISYDPNNPQPLTSRV
VTLWYRPPELLLGATEYGAAVDMWSTGCIVAELFTGKPIMPGRTEVEQIHKIFKLCGSPS
ENYYKKSKVPETAMFKPQQQYRRCVTETFKDLPPSAVLLIDSLLSLEPEVRGTAASALQS
DFFRTKPFACDPSSLPKLPPSKEYDIRLRQEEARRQRNAALGRQGAESIKPGNENHAASR
AIDIAAEVKQPTHNTSKSTCEKFNTEDSVPGFRVEPRALPTSMQVPECGSTWNNTGGYAD
HRSVLGRVYSSVRVARKKGSSNSNIPQYDAADLRNDIEITDHNQQVDRPVSSQKKEQQED
HGRKHKRIHYSGPLMPPGGNIEDMLKEHERHIQEAVRKARLSKGSR"; # Your actual protein string

#my $current_md5 = md5_hex(uc($sequence));
#print "Current Sequence MD5: $current_md5\n";
#3af50c6abad1611232170c88581f8943
#md5 is 3af50c6abad1611232170c88581f8943
#translation start is 472  
#aa change is G -> E, at pos 472
#sequence with respect to the transcript: A
#sequence with respect to the variation feature: A
#offset: 18849, matrix len = 3883
#-------------------- WARNING ----------------------
#MSG: Offset outside of prediction matrix for position 472 and amino acid E?
#FILE: EnsEMBL/Variation/ProteinFunctionPredictionMatrix.pm LINE: 718
#CALLED BY: EnsEMBL/Variation/ProteinFunctionPredictionMatrix.pm  LINE: 349
#Date (localtime)    = Fri Mar 13 11:21:58 2026
#Ensembl API version = 108
#---------------------------------------------------
#consequence SO terms: missense_variant
#amino acid change: G/E
#resulting codon: GAA
#reference codon: GGA
#Use of uninitialized value $sift_pred in concatenation (.) or string at gramene-live/scripts/load-scripts/restore_sift_scores.pl line 155.
#SIFT prediction:
#Use of uninitialized value $sift_score in concatenation (.) or string at gramene-live/scripts/load-scripts/restore_sift_scores.pl line 156.
#SIFT score:

perl  gramene-live/scripts/load-scripts/test_fetch_sift_from_pfpm.pl -registry ~/data/genomes/Sbicolor/SorghumBase/Release9/ensembl.registry -s sorghum_bicolor
prot for SORBI_3008G072400.1 is 
MVICCQLQRLARKTAKIGLIAINSFKKKALHLMGCLCSKGAKDDANATSGRRTPSRKSDSAADAVSNNGGTAVLNAKAKEKLSGGEKVAVALDARISSGNNAELKGLSGEHVVAGWPSWLINVAPKAVEGWLPRRADSFEKLAKIGQGTYSIVYKARDLESGKIVALKKVRFVNMDPESVRFMAREIHILRRLDHPNVIKLEGIVTSRVSQNLYLVFEYMEHDLAGLVATPGLKLTEPQIKCFVQQLLHGLDHCHKNGVLHRDIKGANLLIDSNGMLKIGDFGLAISYDPNNPQPLTSRVVTLWYRPPELLLGATEYGAAVDMWSTGCIVAELFTGKPIMPGRTEVEQIHKIFKLCGSPSENYYKKSKVPETAMFKPQQQYRRCVTETFKDLPPSAVLLIDSLLSLEPEVRGTAASALQSDFFRTKPFACDPSSLPKLPPSKEYDIRLRQEEARRQRNAALGRQGAESIKPGNENHAASRAIDIAAEVKQPTHNTSKSTCEKFNTEDSVPGFRVEPRALPTSMQVPECGSTWNNTGGYADHRSVLGRVYSSVRVARKKGSSNSNIPQYDAADLRNDIEITDHNQQVDRPVSSQKKEQQEDHGRKHKRIHYSGPLMPPGGNIEDMLKEHERHIQEAVRKARLSKGSR
Current Sequence MD5: 3af50c6abad1611232170c88581f8943
offset: 45, matrix len = 3883
pfs: 0xc00a => deleterious - low confidence (0.01)
A mutation to 'C' at position 2 is predicted to be deleterious - low confidence, score 0.01
Use of uninitialized value in concatenation (.) or string at gramene-live/scripts/load-scripts/test_fetch_sift_from_pfpm.pl line 127.
Matrix Peptide Length: 
Matrix Object MD5:    3af50c6abad1611232170c88581f8943
Querying Position: 472
offset: 18849, matrix len = 3883

-------------------- WARNING ----------------------
MSG: Offset outside of prediction matrix for position 472 and amino acid E?
FILE: EnsEMBL/Variation/ProteinFunctionPredictionMatrix.pm LINE: 718
CALLED BY: EnsEMBL/Variation/ProteinFunctionPredictionMatrix.pm  LINE: 349
Date (localtime)    = Mon Mar 16 15:39:01 2026
Ensembl API version = 108
---------------------------------------------------
Use of uninitialized value $prediction in concatenation (.) or string at gramene-live/scripts/load-scripts/test_fetch_sift_from_pfpm.pl line 132.
Use of uninitialized value $score in concatenation (.) or string at gramene-live/scripts/load-scripts/test_fetch_sift_from_pfpm.pl line 132.
A mutation to 'E' at position 472 is predicted to be , score


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
				$_->pep_allele_string
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

