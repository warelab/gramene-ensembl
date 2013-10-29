#!/bin/env/perl -w

=pod

=head1 NAME

perl dump_interpro_annotation_for_trees.pl - Dumps interpro annotation summary for genes of each protein tree

=head1 SYNOPSIS

  perl dump_interpro_annotation_for_trees.pl [options] > resultfile.dat

Options:

  -h|--help             Show brief help and exit.
  -m|--man              Show detailed help
  -e|--ensembl_registry Path to Ensembl registry file.

=head1 OPTIONS

B<-h|--help>
  Print a brief help message and exits.

B<-m|--man>
  Print man page and exit

B<-e|-ensembl_registry>
  Use this Ensembl registry file for database connection info.
  Default is <ENSEMBLHOME>/conf/ensembl.registry

=head1 DESCRIPTION

Script that iterates through each protein tree in the compara database
contained in the registry, and examines the member genes for conserved
interpro domains.

The outputted fields are
  print join( "\t",
              'protein_tree_root_node_id',
              'protein_tree_root_taxon_name',
              'protein_tree_leaf_count',
              ( map{s/\s+/_/g;$_} sort keys %all_species ) );

B<The Ensembl Registry>

  The database connection details for both Ensembl species and compara
  databases must be configured in an ensembl registry file.

  An example registry file would be;
  ---
  use Bio::EnsEMBL::DBSQL::DBAdaptor; 
  use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
  use Bio::EnsEMBL::Registry;
  # The Ensembl compara database
  Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new
  ( '-species' => 'compara',
    '-group'   => 'compara',
    '-dbname'  => 'ensembl_compara_48_28', );
  ---

TODO: Complete this section
                   
Maintained by Will Spooner <whs@ebi.ac.uk>

=cut


use strict;
use warnings;
use Data::Dumper qw(Dumper);
use File::Basename;
use FindBin qw( $Bin );
use Pod::Usage;
use Getopt::Long;
use IO::File;

use vars qw( $BASEDIR );
BEGIN{
  # Set the perl libraries 
  $BASEDIR = dirname( dirname( dirname($Bin) ) );
  -d $BASEDIR.'/ensembl-live'
    || die( "\n[*DIE] Need $BASEDIR/ensembl-live symlinked to ensembl\n" );
  unshift @INC, $BASEDIR.'/ensembl-live/ensembl/modules';
  unshift @INC, $BASEDIR.'/ensembl-live/ensembl-compara/modules';
}

use Bio::EnsEMBL::Registry;

our $ENS_DBA;
our $CMP_DBA;
BEGIN{  #Argument Processing
  my $help=0;
  my $man=0;
  my( $species, $old_species, $reg, $logic_name, $no_insert );
  GetOptions
      ( 
        "help|?"             => \$help,
        "man"                => \$man,
        "ensembl_registry=s" => \$reg,
        )
      or pod2usage(2);
  pod2usage(-verbose => 2) if $man;
  pod2usage(1) if $help;

  $reg        ||= $BASEDIR.'/conf/ensembl.registry';
  map{
    -e $_ || pod2usage( "\n[*DIE]File $_ does not exist\n" );
    -r $_ || pod2usage( "\n[*DIE]Cannot read $_\n" );
    -f $_ || pod2usage( "\n[*DIE]File $_ is not plain-text\n" );
    -s $_ || pod2usage( "\n[*DIE]File $_ is empty\n" );
  } $reg;
  
  # Load the ensembl file
  Bio::EnsEMBL::Registry->load_all( $reg );
  $CMP_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( 'compara','compara')
      || pod2usage("\n[*DIE]No compara DB set in $reg\n" );
}


MAIN:{

  # Collect the species in the DB
  my $genomedb_adaptor  = $CMP_DBA->get_adaptor('GenomeDB');
  my %species;
  my %interpros;
  foreach my $genome_db( @{$genomedb_adaptor->fetch_all} ){
    $species{$genome_db->dbID} = $genome_db->name;
  }
  
  # Create the adaptors that we will be using to fetch data from the database
  my $tree_adaptor  = $CMP_DBA->get_adaptor('ProteinTree');

  # Some counters
  my( $num_trees ) = (0);

  # Degug
  my $db = $CMP_DBA->dbc->dbname;
  warn( "[INFO] Collecting all trees from $db... \n" );

  # Get a list of all protein trees, and loop through them
  my ($root) = @{$tree_adaptor->fetch_all_roots};
  my $tree_arrayref = $root->children;
  my $total_trees = scalar @$tree_arrayref;

  warn( "[INFO] Collected trees $total_trees; processing...\n" );

  # Variables in which we store data for printing
  my %all_species;
  my %all_taxons;
  my %all_interpros;
  my %tree_species_counts;
  my %tree_interpro_counts;


  # loop through each tree
  while( my $tree = shift @$tree_arrayref ){
    $num_trees ++; #COUNT
    my $tree_id = $tree->node_id;
    my $tree_taxon = $tree->get_tagvalue('taxon_name') || '';

    $tree_species_counts{$tree_id}->{TREE_TAXON} = $tree_taxon;
    $tree_species_counts{$tree_id}->{LEAF_COUNT} = 0;
    $tree_species_counts{$tree_id}->{IP_GENE_COUNT} = 0;
    $tree_species_counts{$tree_id}->{IP_DESC} = {};
    $tree_species_counts{$tree_id}->{IP_COUNT} = {};
    $tree_species_counts{$tree_id}->{IP} = {};
    $tree_species_counts{$tree_id}->{SPECIES} = {};

    my $node_arrayref = $tree->get_all_leaves;

    foreach my $node( @$node_arrayref ){ # Loop through each node
      my $species = $species{$node->genome_db_id};
      $tree_species_counts{$tree_id}->{LEAF_COUNT} ++;
      $tree_species_counts{$tree_id}->{SPECIES}{$species} ++;
      $all_species{$species} ++;
      
      my $translation = $node->get_Translation;
      unless( $translation ){
        warn( sprintf("[WARN] Cannot get %s translation %s",
                          $species,$node->stable_id) );
        next;
      }
      my $gene_has_ip;
      foreach my $pf( @{$translation->get_all_ProteinFeatures} ){
        my $ip = $pf->interpro_ac || next;
        $tree_species_counts{$tree_id}{IP_DESC}{$ip} = $pf->idesc;
        $tree_species_counts{$tree_id}{IP}{$ip} ||= {};
        $tree_species_counts{$tree_id}{IP}{$ip}{$species} ++;
        $tree_species_counts{$tree_id}{IP_GENE_COUNT} ++ unless $gene_has_ip;
        $tree_species_counts{$tree_id}{IP_COUNT}{$ip} ++ unless $gene_has_ip->{$ip};
        $gene_has_ip->{$ip} ++;
      }
    }
    $tree->release_tree;

    # Update progress each 1000 genes
    if( $num_trees/1000-int($num_trees/1000) == 0 ){
      warn( "[INFO] processed $num_trees of $total_trees\n" );
    }   
#    if( $num_trees == 10 ){
#      warn Dumper(\%tree_species_counts);
#      last;
#    }
  }

  # print the header
  print join( "\t",
              'root_node_id',
              'root_taxon_name',
              'species_count',
              'gene_count',
              'gene_with_ip_count',
              'most_common_ip',
              'most_common_ip_desc',
              'most_common_ip_gene_count',
              'most_common_ip_species_count',
              'shared_AtOs_ip' );
              
  print "\n";

  # Print the results
  foreach my $tree_id( keys %tree_species_counts ){ 
    #warn Dumper( $tree_species_counts{$tree_id} );
    my $rec = $tree_species_counts{$tree_id};
    my( $most_common_ip ) = ( sort{$rec->{IP_COUNT}->{$b}
                                   <=> $rec->{IP_COUNT}->{$a} } 
                              keys %{$rec->{IP_COUNT}} );
    my $shared_atos = -1; # Assume no Os or At in tree
    if( $rec->{SPECIES}{'Oryza sativa'} 
        and $rec->{SPECIES}{'Arabidopsis thaliana'} ){
      # Have both At and Os in tree
      $shared_atos = 0;
      foreach my $ip( keys %{$rec->{IP}} ){
        if( $rec->{IP}{$ip}{'Oryza sativa'}
            and $rec->{IP}{$ip}{'Arabidopsis thaliana'} ){
          # We have shared InterPro
          $shared_atos ++;
        }
      }
    }
    print join( "\t",
                $tree_id,
                $rec->{TREE_TAXON},
                scalar(keys %{$rec->{SPECIES}}),
                $rec->{LEAF_COUNT},
                $rec->{IP_GENE_COUNT},
                $most_common_ip || '',
                $most_common_ip ? $rec->{IP_DESC}->{$most_common_ip} : 0,
                $most_common_ip ? $rec->{IP_COUNT}->{$most_common_ip} : '',
                $most_common_ip 
                ? scalar( keys %{$rec->{IP}->{$most_common_ip}}) : 0,
                $shared_atos );
    print "\n";
  }

}
exit;

__END__

