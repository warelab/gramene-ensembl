#!/usr/bin/env perl
=head1 NAME

delete_genes_by_logicname.pl - delete all the genes and their transcript and associated attribs.

	

=cut

BEGIN {
    $ENV{'GrameneDir'} = '/usr/local/gramene/'; 
    $ENV{'GrameneEnsemblDir'} = '/usr/local/ensembl-live/'; 
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


=head1 SYNOPSIS

dump_fasta.pl  [options] 
 
 Options:
    --help		help message
    --man		full documentation
    --registry          the registry file for database connections
    --speceis 		which species to dump
    --exclude		genes to exclude
    --bylogicname	make a fasta file for each analysis logic_name
    --analysispgm	this(these) analysispgm(s) only
    --exclude-analysispgm	skip this analysispgm
    --exclude-clone	skip this clone (when doing all genes)


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

supply the species name whose genes/transcripts are to be deleted


=item B<--help> 

print a help message and exit

=item B<--man> 

print documentation and exit

    
=back

=head1 ARGUMENTS

Gene Stable ids - only dump transcripts of these genes
None=All.

=cut

my ($species, $registry, @types);
my (%exclude_gene,%exclude_analysispgm,%analysispgm,%exclude_clone,%bylogicname);

{  #Argument Processing
  my $help=0;
  my $man=0;
  my @exclude_gene=();
  my @exclude_analysispgm=();
  my @analysispgm=();
  my @exclude_clone=();
  my @bylogicname=();
  GetOptions( "help|?"=>\$help,"man"=>\$man
	      ,"exclude=s"=>\@exclude_gene
	      ,"bylogicname=s"=>\@bylogicname
	      ,"exclude-analysispgm=s"=>\@exclude_analysispgm
	      ,"analysispgm=s"=>\@analysispgm
	      ,"exclude-clone=s"=>\@exclude_clone
	      ,"species=s"=>\$species
	      ,"registry=s"=>\$registry

	    )
    or pod2usage(2);
  pod2usage(-verbose => 2) if $man;
  pod2usage(1) if $help;
  #pod2usage(2) if $margin<0;
  %exclude_gene= map { $_,1 }  map { split /,/ } @exclude_gene;
  %exclude_analysispgm= map { $_,1 } map { split /,/ } @exclude_analysispgm;
  %analysispgm= map { $_,1 }  map { split /,/ } @analysispgm;
  %bylogicname= map { $_,1 }  map { split /,/ } @bylogicname;
  %exclude_clone= map { $_,1 } map { split /,/ } @exclude_clone;
}


# Load the ensembl file
Bio::EnsEMBL::Registry->load_all( $registry );
my $ENS_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'core' );
$ENS_DBA || pod2usage( "\nNo core DB for $species set in $registry\n" );

my $gene_adaptor=$ENS_DBA->get_GeneAdaptor;
my $slice_adaptor = $ENS_DBA->get_SliceAdaptor; 
#my $assembly_type=$ENS_DBA->get_MetaContainer()->get_default_assembly();

my @genes=@ARGV;
@genes or  @genes=@{$gene_adaptor->list_dbIDs()};

my %count;

my (%seqio,%seqio5,%seqio3);


foreach my $geneid (@genes) {
  #print "! $geneid\n";
  
  $count{total_genes}++;
  my $gene;
  if(@ARGV) {    # Stable ids GRMGnnnnnnn
    eval { $gene=$gene_adaptor->fetch_by_stable_id($geneid) };
    print STDERR "$@\n" and next if $@;
  } else {	#internal ids
    eval { $gene= $gene_adaptor->fetch_by_dbID($geneid); };
    print STDERR "gene_id $geneid:\n$@\n" and next if $@;
    # fails e.g. if gene has no exons
  }
  next unless $gene;
  
  next if %exclude_gene and $exclude_gene{$gene->stable_id}
    or %analysispgm and !$analysispgm{$gene->analysis->program}
    or %bylogicname and !$bylogicname{$gene->analysis->logic_name} 
      or %exclude_analysispgm and $exclude_analysispgm{$gene->analysis->program}
	;
  
  if (%exclude_clone) {
    my $slice;
    eval {	#have some genes on phase1 clones -- no assembly, no slice
      $slice=$slice_adaptor->fetch_by_gene_stable_id($gene->stable_id,0);
    };
    print STDERR "$@\n" and next if $@;
    foreach my $tile (@{$slice->get_tiling_path}) {
      next GENE if $exclude_clone{$tile->component_Seq->clone->embl_id};
    }
    
  }
  
  $count{qualified_genes}++;
  $gene_adaptor->remove($gene);
  
  
}

for my $k (sort keys %count){
    print "$k = $count{$k}\n";
}




