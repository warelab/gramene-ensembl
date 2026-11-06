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

use Bio::EnsEMBL::Variation::TranscriptVariationAllele;
#$tva->sift_prediction
#$tva->sift_score

#use List::Compare;

=head1 SYNOPSIS

dump_transcripts.pl  [options] 
 
 Options:
    --help		help message
    --man		full documentation
    --registry          the registry file for database connections
    --speceis 		which species to dump
    --debug		debug
    --format		two types, 1 (default) or 2 (dump actual src:sample_name)
    --chr		dump by chr
 
=head1 OPTIONS

=over 4

=item B<--registry>

The registry file for ensembl databases

=item B<--species> 

supply the species name whose transcripts are to be dumped

=item B<--chr>

only dump varaiation for this chr. 

=item B<--format>

reporting format
gt_sum: print out the count of samples for each genotype
	Transcript(strand)	Allele(transcript_allele)	ConseqType	PositionInTranscript	PositionInCDS	PositionInProt	AminoAcid	Codons
gt_lines print out the actaul sample for each genotype in the form of GT1:source:sample|source:sample|...|source:sample,GT2:source:sample|source:sample|...|source:sample

=item B<--help> 

print a help message and exit

=item B<--man> 

print documentation and exit

    
=back


=cut

our @LOF = qw(
Splice_Acceptor_Variant
Splice_donor_Variant
Stop_Gained
Frameshift_Variant
Stop_Gained
Start_Lost
Missense_Variant
);

our $LOF_EXP = lc (join '|', sort @LOF);
 


# Osmh63.01G000010_02.1 as example
# format 1
#Transcript (strand)	Allele (transcript allele)	Consequence Type	Position in transcript	Position in CDS	Position in protein	Amino acid	Codons
#format 2
#Variant ID	Chr: bp	vf_allele	Alleles	Conseq. Type	Transcript

our ($species, $registry, $format, $debug,  $CHR);

{  #Argument Processing
  my $help=0;
  my $man=0;
  $format = 1;

  GetOptions( "help|?"=>\$help,"man"=>\$man
	      ,"species=s"=>\$species
	      ,"registry=s"=>\$registry
	      ,"format=i"=>\$format
	      ,"debug"=>\$debug
	      ,"chr=s" => \$CHR
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
our $sgta =  $VAR_DBA->get_adaptor('samplegenotype');
# my $sgta = $registry->get_adaptor('human', 'variation', 'samplegenotype');

my %count;
my $slices = $slice_adaptor->fetch_all('toplevel');
our $cnt = 0;
	
for ( @$slices ){ 
	my @transcripts;
	my @print_out;
	my $seq_name = $_->seq_region_name ;
#	warn ("CHR=$CHR and seq_name=$seq_name\n");
	next if ( (defined $CHR) and (lc ($CHR) ne lc ($seq_name)) );
	warn("Process chr $seq_name\n"); 
	my $outfile = "TrptVar.$seq_name.report";
	push @transcripts, 

	map{
           		my $id = $_->stable_id;
             		#print "! $id\n" if $debug;
              		push @print_out, map{join "\t", @{$_}}
				 @{get_transcript_variations($_)};
              		$count{trpt}{$seq_name}++;
			map{ print "$_\n";} @print_out and exit if ($debug && $cnt > 5);

	}@{$_->get_all_Transcripts};

	warn("About to print to $outfile\n");	
	$count{trpt_var}{$seq_name} += scalar @print_out; 
	open my $fh, ">$outfile" or die "Cannot open $outfile to write";
       	for ( @print_out){
              		print $fh "$_\n";
       	}
       	close $fh;
	last if $debug ;
}




for my $k (sort keys %count){
	my $total = 0;
      	for my $chr( keys %{$count{$k}}){
		$total += $count{$k}{$chr};	
            	print STDERR "$k $chr  = $count{$k}{$chr}\n";
      	}
	print STDERR "$k total = $total\n";
}

sub get_transcript_variations {
  	my $transcript   = shift;

	my @tvs = @{ $tv_adaptor->fetch_all_by_Transcripts([$transcript]) };
	my $tvs_cnt = scalar @tvs;
	print "Get $tvs_cnt transcript variation pairs for " . $transcript->stable_id ."\n" if $debug;

	my @tvs_report;
	my $lof_report;
	for my $tv ( @tvs ){
                my $vf = $tv->variation_feature();
       
		my $variation_name = $vf->variation_name;
		next unless $variation_name =~ /^rs/i;

		my $transcript_stable_id = $tv->transcript_stable_id;
		my $gene_id = get_gene_id($transcript_stable_id); 
		my $conseq = $tv->display_consequence; 
		my $chr    = $vf->seq_region_name;
		my $seq_region_start = $vf->seq_region_start;
		my $seq_region_end   = $vf->seq_region_end;
		my $seq_region_strand   = $vf->seq_region_strand;
		my @alt_alleles = @{$vf->alt_alleles};
		my $allele_string = $vf->allele_string;
		my $alt_alleles_str = join "|", @alt_alleles;
		my $cdna_start = $tv->cdna_start;
		my $cds_start  = $tv->cds_start;
                my $translation_start = $tv->translation_start;
                my $pep_allele_string = $tv->pep_allele_string;
                my $codons            = $tv->codons;
		
		my $variation = $vf->variation;
		my $source    = $variation->source_name;
                #my @alleles = @{$variation->get_all_Alleles()};   
                
		#my $lc = List::Compare->new(\@alt_alleles, \@alleles);
		#my @ref_allele = $lc->get_complement;


		my @common_fields = (
			$gene_id,
			$transcript_stable_id,
			$variation_name,		
			$source,
			"$chr($seq_region_strand):$seq_region_start-$seq_region_end",
			$allele_string,
			$alt_alleles_str,
			$conseq,
			$cdna_start,
			$cds_start,
			$translation_start,
			$pep_allele_string,
			$codons,
		);


#Osmh63.01G000010_02     Chr01_11742_C_G 1(1):11742-11742        C/G     ARRAY(0x4343f60)        ARRAY(0x4340df0)        missense_variant 790     415     139     P/A     Cca/Gca       
	        
		if( $conseq =~ /^($LOF_EXP)$/i){
			#find samples with deleterious alle

			my %genotype2samples;
			my $sample_genotypes = $variation->get_all_SampleGenotypes();
  			foreach my $sgt (@$sample_genotypes) {
    				my $sgt_genotype = $sgt->genotype_string;
    				my $sgt_sample_name = $sgt->sample()->name;
    				print "$variation_name\t$sgt_genotype\t$sgt_sample_name\n" if $debug;
	
				push @{$genotype2samples{$sgt_genotype}},$sgt_sample_name;
			}
			
			#calculate sift score if possible
			# get TranscriptVariationAlleles for sift score
                        my $tva = $tv->get_all_TranscriptVariationAlleles();
                        my $sift_scores = join ',', (map{ $_->sift_score } @$tva);


			#get 500bp flanking seq
			#my ($upstream_flank_start, $upstream_flank_end) = ($seq_region_start-250, $seq_region_start-1);
			#my ($downstream_flank_start, $downstream_flank_end) = ($seq_region_end+1, $seq_region_end+250);
			#
			#my ($upstream_seq, $downstream_seq);
			#
			#if( $seq_region_strand == 1 ){
			#	$upstream_seq = $slice_adaptor->fetch_by_location("$chr:$upstream_flank_start-$upstream_flank_end", 'toplevel')->seq;
			#	$downstream_seq = $slice_adaptor->fetch_by_location("$chr:$downstream_flank_start-$downstream_flank_end", 'toplevel')->seq; 				
			#}else{
			#	 warn("[WARN] There are variations on - strand\n");
			#
			#}


			my @sorted_gt_by_samplecnt = sort { scalar @{$genotype2samples{$a}} <=> scalar @{$genotype2samples{$b}} }
                               (keys %genotype2samples);
			my @gt_samples; 
			if( $format == 1){
				warn("DEBUG format is $format\n") && $cnt++ if $debug;
				push @gt_samples, map { $_ . ':' . scalar @{$genotype2samples{$_}} } @sorted_gt_by_samplecnt;
			
			}elsif( $format == 2 ){
				push @gt_samples, map { $_ . ':' . join '|', @{$genotype2samples{$_}} } @sorted_gt_by_samplecnt;
			}else{ warn "format $format does not exist";}

			push @tvs_report,[@common_fields, $sift_scores, 
						(join ',', @gt_samples), 
					];


		}else{
#print "Not LOF\n" if $debug;
#			push @tvs_report,[@common_fields];
        	}

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
				get_gene_id($_->transcript_stable_id).':'.$_->transcript_stable_id,
				$_->pep_allele_string, 
				$_->translation_start,
				]
				} @{$tvs};
		push @snp_report, @asnp_tvs;				
	}

     	return \@snp_report;
}


sub get_gene_id {

	my $trpd_id= shift;
	my $gene_stable_id = $trpt_adaptor->fetch_by_stable_id($trpd_id)
				->get_Gene()->stable_id;
	
	return $gene_stable_id; 


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

