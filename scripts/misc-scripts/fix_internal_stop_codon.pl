#!/usr/local/bin/perl 

=head1 NAME

fix_internal_start_codon.pl - Some translation start before Methonine, Find all the protein coding genes, check if they all start from M, if not, find the 1st M and set start codon there.
	

=cut

BEGIN {
    $ENV{'GrameneDir'} = '/usr/local/gramene/'; 
    $ENV{'GrameneEnsemblDir'} = '/usr/local/ensembl-live/'; 
}

# The first shall be last...
use lib map { $ENV{'GrameneDir'}."/$_" } qw ( lib/perl );

use lib map { $ENV{'GrameneEnsemblDir'}."/$_" } 
        qw ( bioperl-live modules ensembl/modules conf
	     ensembl-external/modules ensembl-draw/modules
	     ensembl-compara/modules );

use strict;
use warnings;
use Scalar::Util qw (reftype);

#use Gramene::Config;
use Bio::SeqIO;
use Bio::EnsEMBL::DBLoader;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::TranscriptMapper;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Exon;

use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

=head1 SYNOPSIS

fix_internal_stop_codon.pl  [options] 
 
 Options:
    --help		help message
    --man		full documentation
    --registry          the registry file for database connections
    --speceis 		which species to dump
    --exclude		genes to exclude
    --bylogicname	only work on genes udner this analysis logic_name
    --nowrite           test only , no change to the database
    --debug             debug mode

=head1 OPTIONS

=over 4


=item B<--exclude>

    Genes to ignore.  

=item B<--bylogicname>

    Only work on Genes under this logicname

=item B<--registry>

    The registry file for ensembl databases

=item B<--species> 

    supply the species name whose transcripts are to be dumped

=item B<--nowrite>

    test run without writing to database

=item B<--help> 

    print a help message and exit

=item B<--man> 

    print documentation and exit

=item B<--debug> 

   print out more debug information

=back

=head1 ARGUMENTS

    Gene Stable ids - only dump transcripts of these genes
    None=All.

=cut

my ($species, $registry);
my (%exclude_gene, $bylogicname, $debug, $nowrite);
my $dump_filename='internalStoptranscript_stable_id';
{  #Argument Processing
  my $help=0;
  my $man=0;
  my @exclude_gene=();

  GetOptions( "help|?"=>\$help,"man"=>\$man
	      ,"exclude=s"=>\@exclude_gene
	      ,"bylogicname=s"=>\$bylogicname
	      ,"species=s"=>\$species
	      ,"registry=s"=>\$registry
	      ,"debug"=>\$debug
	      ,"nowrite"=>\$nowrite
	    )
    or pod2usage(2);
  pod2usage(-verbose => 2) if $man;
  pod2usage(1) if $help;
  %exclude_gene= map { $_,1 }  map { split /,/ } @exclude_gene;
 
}


# Load the ensembl file
Bio::EnsEMBL::Registry->load_all( $registry );
my $ENS_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'core' );
$ENS_DBA || pod2usage( "\nNo core DB for $species set in $registry\n" );

my $gene_adaptor=$ENS_DBA->get_GeneAdaptor;
my $aa = $ENS_DBA->get_AttributeAdaptor();

print "DBbase connected is ", $ENS_DBA->dbname, "\n" if $debug;

my @genes = map{ $gene_adaptor->fetch_by_stable_id($_) } @ARGV;

#my @genes = ($gene_adaptor->fetch_by_stable_id('Opunc01g00010'));
@genes or  
    @genes = $bylogicname ? @{$gene_adaptor->fetch_all()}:
    @{$gene_adaptor->fetch_all_by_logic_name()};

my %count;
my %transcript_with_internal_stop;

foreach my $gene(@genes) {
#  print "geneid = ", $gene->stable_id, "\n";
  
  $count{total_genes}++;
  
  next if %exclude_gene and $exclude_gene{$gene->stable_id};
  
  $count{qualified_genes}++;
  
  my @transcripts;
  @transcripts = @{$gene->get_all_Transcripts};

  for my $transcript (@transcripts) {

    my $id = $transcript->dbID;
    my $stableid = $transcript->stable_id;
    my $slice_name = $transcript->slice->name;
    my $biotype = $transcript->biotype;
    my $logic_name = $transcript->analysis->logic_name;
    my $strand = $transcript->strand;
    my $comp_id = join "|", ($id, $stableid, $strand, $slice_name, $logic_name, $biotype);

    my $translation = $transcript->translation();
    my $aa_seq;

    eval{$aa_seq= $translation->seq()};

    if($@){
	print "error translation seq for $comp_id, $@\n";
	next;
    }
    
    if($aa_seq){
	$count{qualified_transcripts}++;
    }else{
	warn("no aa translation for $comp_id\n");
	next;
    }
    
    my $idx = index($aa_seq,'*',0);
    while($idx!=-1) {
	$transcript_with_internal_stop{$stableid} = 1;
	my $pidx = $idx+1;
	$aa->store_on_Translation($translation,
				  [
				   Bio::EnsEMBL::Attribute->new(
				       -CODE        => "amino_acid_sub",
				       -NAME        => "Amino acid substitution",
				       -DESCRIPTION => "Some translations have been manually curated for amino acid substitiutions. For example a stop codon may be changed to an amino acid in order to prevent premature truncation, or one amino acid can be substituted for another.",
				       -VALUE => "$pidx $pidx X")]) unless $nowrite;
	$idx = index($translation->seq(),'*',$idx+1);
    }
  }
  
}

$count{transcript_with_internal_stop} = scalar keys %transcript_with_internal_stop;

for my $k (sort keys %count){
  print "$k = $count{$k}\n";
}

open my $fh, '>', $dump_filename or die "cannot open $dump_filename to write";

print $fh join "\n", keys %transcript_with_internal_stop;


  
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

