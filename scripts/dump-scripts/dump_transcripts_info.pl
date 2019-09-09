#!/usr/local/bin/perl 

=head1 NAME


dump_transcripts_info.pl - dump the info for transcripts 

The format of output is

Bin #	Name	Chr #	Strand	Transcription Start	Transcription End	CDS Start 	CDS End	# Exons	Exon start positions	Exon End Positions	ID	Name 2	cdsStartStat	cdsEndStat	exonFrames

	Sb01g000200	chr1	+	2164	2829	2164	2829	1	2164,	2829,	0	Sb01g000200			0,
	
1643	NM_016459	chr5	-	138751155	138753504	138751352	138753444	4	138751155,138751608,138752048,138753267,	138751509,138751719,138752173,138753504,	0	MGC29506	cmpl	cmpl	2,2,0,0,


=cut


BEGIN {
  $ENV{'GrameneDir'} = '/usr/local/gramene/'; 
  $ENV{'GrameneEnsemblDir'} = '/usr/local/ensembl-53/'; 
}

# The first shall be last...
use lib map { $ENV{'GrameneDir'}."/$_" } qw ( lib/perl );
#use lib map { $ENV{'GrameneEnsemblDir'}."/ensembl-live/$_" } 
#        qw ( bioperl-live modules ensembl/modules ensembl-external/modules
#             ensembl-draw/modules ensembl-compara/modules );
use lib map { $ENV{'GrameneEnsemblDir'}."/$_" } 
        qw ( bioperl-live modules ensembl/modules conf
	     ensembl-external/modules ensembl-draw/modules
	     ensembl-compara/modules );

use strict;
use warnings;


#use Gramene::Config;
use Bio::SeqIO;
use Bio::EnsEMBL::DBLoader;
use Bio::EnsEMBL::Registry;
#use Bio::EnsEMBL::Gene;
use Getopt::Long;
use Pod::Usage;

use Data::Dump qw(dump);

=head1 SYNOPSIS

dump_transcripts_info.pl  [options] 
 
 Options:
    --help		help message
    --man		full documentation
    --registry          the registry file for database connections
    --speceis 		which species to dump
    --bylogicname	make a fasta file for each analysis logic_name "GeneModel_JGI"


=head1 OPTIONS

=over 4


=item B<--bylogicname>

    gene set logic name such as   "GeneModel_JGI"

    
=item B<--registry>

The registry file for ensembl databases

=item B<--species> 

supply the species name whose transcripts are to be dumped


=item B<--help> 

print a help message and exit

=item B<--man> 

print documentation and exit

    
=back

=head1 ARGUMENTS

Gene Stable ids - only dump transcripts of these genes
None=All.

=cut

my ($species, $registry);
my ($bylogicname);

{  #Argument Processing
  my $help=0;
  my $man=0;

  GetOptions( "help|?"     =>\$help,
	      "man"        =>\$man,
	      "bylogicname"=>\$bylogicname,
	      "species=s"  =>\$species,
	      "registry=s"=>\$registry,
	    )
    or pod2usage(2);
  pod2usage(-verbose => 2) if $man;
  pod2usage(1) if $help;

}


# Load the ensembl file
Bio::EnsEMBL::Registry->load_all( $registry );
my $ENS_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'core' );
$ENS_DBA || pod2usage( "\nNo core DB for $species set in $registry\n" );

my $gene_adaptor=$ENS_DBA->get_GeneAdaptor;
my $slice_adaptor = $ENS_DBA->get_SliceAdaptor; 

my %count;

foreach my $gene ( @{$gene_adaptor->fetch_all_by_logic_name( $bylogicname )} ) {
  #print "! $geneid\n";
  
  $count{total_genes}++;
    
  my ($working_set_attrib) = @{$gene->get_all_Attributes('working-set')};
  #print "gene attrib is ", $working_set_attrib->code, "\n";
  next if ($species =~ /zea|mays|maize/i && !$working_set_attrib);
  
  $count{qualified_genes}++;  
  
  foreach my $trans ( @{$gene->get_all_Transcripts()} ) {
  
    #print join "\t", ($trans->stable_id, $trans->spliced_seq, "\n");
    
    my $stable_id          = $trans->stable_id;
    my $seq_region_name    = $trans->seq_region_name;
    my $genomic_start      = $trans->start;
    my $genomic_end        = $trans->end;
    my $coding_start       = $trans->coding_region_start;
    my $coding_end         = $trans->coding_region_end;
    my $strand             = $trans->strand;
    my $exons              = $trans->get_all_Exons;
    

    my @exon_infos         = map  { #print $_->start(), "\n"; 
				    [$_->start(), $_->end(), $_->phase() ] } 
                             sort { $a->start <=> $b->start } @{$exons};     

    #print dump @exon_infos;
    my $cdsStartStat       = $coding_start ? 'cmpl' : 'incmpl';
    my $cdsEndStat         = $coding_end   ? 'cmpl' : 'incmpl';

    my $out = join "\t", (
			  $stable_id,
			  $seq_region_name,
			  $strand,
			  $genomic_start,
			  $genomic_end,
			  $coding_start,
			  $coding_end,
			  $cdsStartStat,
			  $cdsEndStat,
			  scalar @exon_infos,
			  (join ',', (map { $_->[0] } @exon_infos)),
			  (join ',', (map { $_->[1] } @exon_infos)),
			  (join ',', (map { $_->[2] } @exon_infos)),
			 ); 
    
    $count{qualified_transcripts}++;
    
    print "$out\n";
  }
  #last;
}



=head1 OUTPUT


=head1 AUTHOR

   Sharon Wei <weix@cshl.edu>

   Gramene Project (www.gramene.org)
   Ware Lab
   Cold Spring Harbor Laboratory

   This is free software and may be distributed on the same terms as Perl itself.

=cut

