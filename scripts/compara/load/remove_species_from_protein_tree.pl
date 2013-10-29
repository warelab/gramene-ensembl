#!/bin/env perl -w

# This script can be used to remove a specific species from the
# protein tree build

use strict;

# Some constants
#our $REGFILE = '/usr/local/gramene/conf/ensembl.registry';
our $REGFILE = '/home/whs/paper/ensembl.paper.registry';
our $SPECIES = 'Oryza_glaberrima';
our $TAXON   = $SPECIES;
$TAXON =~ s/_/ /g;
our $COMPARA = 'compara';

BEGIN{

  our $BASEDIR = '/usr/local/ensembl-live';

  # Set the perl libraries 
  -d $BASEDIR
    || die( "\n[*DIE] Need a path to the ensembl distribution\n" );
  unshift @INC, $BASEDIR.'/ensembl-live/ensembl/modules';
  unshift @INC, $BASEDIR.'/ensembl-live/ensembl-compara/modules';
}
use Bio::EnsEMBL::Registry;

unless( -f $REGFILE ){
  die( "\n[*DIE] Cannot find ensembl registry $REGFILE\n" );
}

Bio::EnsEMBL::Registry->load_all( $REGFILE );
our $CMP_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $COMPARA,'compara')
    || die("\n[*DIE] No compara DB for $COMPARA set in $REGFILE\n" );

my $genomedb_adaptor = $CMP_DBA->get_adaptor('GenomeDB');
my $mlss_adaptor     = $CMP_DBA->get_adaptor('MethodLinkSpeciesSet');
my $homology_adaptor = $CMP_DBA->get_adaptor('Homology');
my $member_adaptor   = $CMP_DBA->get_adaptor('Member');
my $tree_adaptor     = $CMP_DBA->get_adaptor('ProteinTree');

my $genome_db = $genomedb_adaptor->fetch_by_registry_name( $SPECIES );
$genome_db || die( "\n[*DIE] No species $SPECIES in compara.\n" );

my $genome_db_id = $genome_db->dbID;

foreach my $root( @{$tree_adaptor->fetch_all_roots||[]} ){
  #my @trees = $tree_adaptor->fetch_node_by_node_id('326082');
  my @trees = @{$root->children};
  foreach my $tree( @trees ){

    #warn( 'Tree: ', $tree->node_id );
    #warn( 'Leaves before: ', scalar(@{$tree->get_all_leaves}) );

    my @leaves = @{$tree->get_all_leaves};
    my @sp_leaves 
        = grep{ !  $_->can("genome_db")
                or $_->genome_db->dbID == $genome_db_id } @leaves;

    my @disavowed;
    if( scalar(@leaves) == scalar(@sp_leaves) ){ # Single species tree
      warn( "$TAXON-only tree ROOT: ", $tree->_root_id," > Delete" );
      $tree_adaptor->delete_node_and_under($tree);
    }    
    else{ # Multi-species tree
      foreach my $leaf (@sp_leaves) {
        warn( '  x> ', $leaf->name," ID:",$leaf->node_id
              ," ROOT:",$leaf->_root_id
              ," > disavowing parent\n" );
        $leaf->disavow_parent;
        $tree = $tree->minimize_tree;
        push @disavowed, $leaf;
      }
    }
    if ($tree->get_child_count == 1) {
      my $child = $tree->children->[0];
      $child->parent->merge_children($child);
      $child->disavow_parent;
    }
    #warn( 'Leaves after: ', scalar(@{$tree->get_all_leaves}), "\n\n" );
      
    if( @disavowed ){    
      $tree_adaptor->sync_tree_leftright_index($tree);
      my $new_node = $tree_adaptor->store($tree);
      $tree_adaptor->delete_nodes_not_in_tree($tree);
    }
    $tree->release_tree;
  }
}

1;

__END__

