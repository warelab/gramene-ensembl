#!/usr/local/bin/perl 

=head1 NAME

cleanup_overkap_fgenesh - when fgenesh genes overlap each other, keep only the longest one
	

=cut


#use

#fetch_all_by_logic_name('foobar');

BEGIN {
    $ENV{'GrameneDir'} = '/usr/local/gramene/'; 
    $ENV{'GrameneEnsemblDir'} = '/usr/local/ensembl-live/'; 
}

# The first shall be last...
use lib map { $ENV{'GrameneDir'}."/$_" } qw { lib/perl };

use lib map { $ENV{'GrameneEnsemblDir'}."/$_" } 
        qw { bioperl-live modules ensembl/modules conf
	     ensembl-external/modules ensembl-draw/modules
	     ensembl-compara/modules };

use strict;
use warnings;
use List::MoreUtils qw{ distinct };


use Bio::SeqIO;
use Bio::EnsEMBL::DBLoader;
use Bio::EnsEMBL::Registry;
use Getopt::Long;
use Pod::Usage;
use Readonly;

Readonly my $remove => 'remove';
Readonly my $keep => 'keep';

=head1 SYNOPSIS

dump_gene_coord.pl  [options] 
 
 Options:
    --help		help message
    --man		full documentation
    --registry          the registry file for database connections
    --speceis 		which species to dump
    --logic_name        logic name of the fgenesh analysis

=head1 OPTIONS

=over 4


=item B<--registry>

The registry file for ensembl databases

=item B<--species> 

supply the species name whose transcripts are to be dumped

=item B<--logic_name> 

 logic name of the fgenesh analysis

=item B<--d> 

 delete genes if set

=item B<--help> 

print a help message and exit

=item B<--man> 

print documentation and exit

    
=back

=head1 ARGUMENTS

Input files with Gene Stable ids 

=cut

my ($species, $registry, $logic_name, $delete);

{  #Argument Processing
  my $help=0;
  my $man=0;

  GetOptions( "help|?"=>\$help,
	      "man"=>\$man,
	      "species=s"=>\$species,
	      "registry=s"=>\$registry,
	      "logic_name=s"=>\$logic_name,
	      "d"=>\$delete,
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

my $cnter=0;
foreach my $gene (@{$gene_adaptor->fetch_all_by_logic_name($logic_name)}) {

    $gene_adaptor->remove($gene);
    $cnter++;
 }

print "deleted $cnter genes";

sub print_genes{

    my $genelist=shift;
    my $fh=shift;

    for my $g (@$genelist){
	print $fh (join "\t", (
				$g->dbID,
				$g->stable_id,
				$g->seq_region_name,
				$g->seq_region_start,
				$g->seq_region_end,
				$g->seq_region_strand,
				));
	print $fh "\n";
    }
}





=head1 AUTHOR

   Sharon Wei <weix@cshl.edu>

   Gramene Project (www.gramene.org)
   Ware Lab
   Cold Spring Harbor Laboratory

   This is free software and may be distributed on the same terms as Perl itself.

=cut

