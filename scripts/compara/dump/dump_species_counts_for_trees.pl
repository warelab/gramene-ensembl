#!/bin/env/perl -w

=pod

=head1 NAME

perl dump_species_counts_for_trees.pl - Dumps species counts for each protein tree

=head1 SYNOPSIS

  perl dump_species_counts_for_trees.pl [options] > resultfile.dat

Options:

  -h|--help             Show brief help and exit.
  -m|--man              Show detailed help
  -e|--ensembl_registry Path to Ensembl registry file.
  -i|--internal_nodes   Dump internal nodes

=head1 OPTIONS

B<-h|--help>
  Print a brief help message and exits.

B<-m|--man>
  Print man page and exit

B<-e|-ensembl_registry>
  Use this Ensembl registry file for database connection info.
  Default is <ENSEMBLHOME>/conf/ensembl.registry

B<-i|-internal_nodes>
  Dump counts for internal nodes of each tree in addition to the 
  leaf nodes.

=head1 DESCRIPTION

Script that iterates through each protein tree in the compara database
contained in the registry, and counts the number of genes for each
species.

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

use vars qw( $BASEDIR $INTERNAL );
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
        "internal_nodes"     => \$INTERNAL,
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
  my %tree_species_counts;
  my %tree_taxon_counts;

  # loop through each tree
  while( my $tree = shift @$tree_arrayref ){
    $num_trees ++; #COUNT
    my $tree_id = $tree->node_id;
    my $tree_taxon = $tree->get_tagvalue('taxon_name') || '';

    $tree_species_counts{$tree_id}->{TREE_TAXON} = $tree_taxon;
    $tree_species_counts{$tree_id}->{LEAF_COUNT} = 0;
    $tree_taxon_counts{$tree_id}->{NODE_COUNT} = 0;

    my $node_arrayref = $INTERNAL # Don't bother with internal nodes unless 
        ? $tree->get_all_nodes    # requested
        : $tree->get_all_leaves;

    foreach my $node( @$node_arrayref ){ # Loop through each node
      if( $node->is_leaf ){
        unless( $node->can('gene_member') ){
          warn sprintf( "$node %s cannot->gene_member", $node->node_id );
          next;
        }
        my $species = $species{$node->gene_member->genome_db_id};
        $tree_species_counts{$tree_id}->{LEAF_COUNT} ++;
        $tree_species_counts{$tree_id}->{$species} ++;
        $all_species{$species} ++;
      }
      else { # Internal node
        my $taxon = $node->get_tagvalue('taxon_name') || '';
        $tree_taxon_counts{$tree_id}->{NODE_COUNT} ++;
        $tree_taxon_counts{$tree_id}->{$taxon} ++;
        $all_taxons{$taxon} ++;
      }
    }
    $tree->release_tree;

    # Update progress each 1000 genes
    if( $num_trees/1000-int($num_trees/1000) == 0 ){
      warn( "[INFO] processed $num_trees of $total_trees\n" );
    }    
  }

  # print the header
  print join( "\t",
              'protein_tree_root_node_id',
              'protein_tree_root_taxon_name',
              'protein_tree_leaf_count',
              ( map{s/\s+/_/g;$_} sort keys %all_species ) );
  if( $INTERNAL ){
    print "\t";
    print join( "\t",
                'protein_tree_internal_node_count',
                ( map{s/\s+/_/g;$_} sort keys %all_taxons ) );
  }
  print "\n";

  # Print the results
  foreach my $tree_id( keys %tree_species_counts ){ 
    #warn Dumper( $tree_species_counts{$tree_id} );
    print join( "\t",
                $tree_id,
                $tree_species_counts{$tree_id}->{TREE_TAXON},
                $tree_species_counts{$tree_id}->{LEAF_COUNT},
                ( map{$tree_species_counts{$tree_id}->{$_}||0} 
                  sort keys %all_species ) );
    if( $INTERNAL ){
      print "\t";
      print join( "\t",
                  $tree_taxon_counts{$tree_id}->{NODE_COUNT},
                  ( map{$tree_taxon_counts{$tree_id}->{$_}||0} 
                    sort keys %all_taxons ) );      
    }
    print "\n";
  }

}
exit;

__END__


    # Get the subnode for the gene tree
    my $subroot;
    if( my $member  
        =  $member_adaptor->fetch_by_source_stable_id('ENSEMBLGENE',$id) ){

      if( my $gene_tree 
          = $tree_adaptor->fetch_by_Member_root_id($member, 0) ){

        $subroot   = $gene_tree->subroot;
        $num_tree_genes                ++; #COUNT
        $tree_nodes{$subroot->node_id} ++; #COUNT
      }
    }

    # Print the data
    print join( "\t",
                $id,
                $gene->coord_system_name,
                $gene->seq_region_name,
                $gene->seq_region_start,
                $gene->seq_region_end,
                $subroot 
                ? ( $subroot->node_id,
                    $subroot->get_tagvalue('taxon_name') || '' )
                : ( '','' ) );
    print "\n";

    # The following lists all node annotations if you want to add extra
    # my @tags = $subroot->get_all_tags;
    # warn join( ", ", @tags );

    # Update progress each 1000 genes
    if( $num_genes/1000-int($num_genes/1000) == 0 ){
      warn( "[INFO] processed $num_genes of $total_genes\n" );
    }
    
  }

  # Debug
  warn sprintf( "[INFO] processed %d coding genes; %d in %d unique trees\n",
                $num_coding, $num_tree_genes, scalar( keys( %tree_nodes ) ) ); 

}

exit;
