#!/bin/env perl -w

# A quick-and-dirty script to identify species specific clusters in
# gene trees.
# First we get a list of all paralogs for the given species, then,
# for each gene in the paralog pair, if we have not seen the gene
# before (genes are unique in all trees):
#  Get the tree node for the gene,
#  Descend to the parent node, and see if all leaves are the same species,
#  If true, then recurse.
#  If false, then go back to the previous node,
#  Dump the list of genes for this sub tree
#
# This script has to be run from brie to get the paths correct

use strict;

# Some constants
our $REGFILE = '/usr/local/gramene/conf/ensembl.registry';
our $SPECIES = 'Zea_mays';
our $COMPARA = 'compara';
our $METHOD_LINK = 'ENSEMBL_PARALOGUES';
our $OUTFILE = "${SPECIES}_specific_clusters.tsv";
our @FIELDS  = qw( tree_node_id
                   tree_leaf_count
                   tree_taxon_name
                   cluster_node_id
                   cluster_leaf_count
                   cluster_taxon_name
                   cluster_genes );

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
our $ENS_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $SPECIES, 'core' )
    || die("\n[*DIE] No core DB for $SPECIES set in $REGFILE\n" );
our $CMP_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $COMPARA,'compara')
    || die("\n[*DIE] No compara DB for $COMPARA set in $REGFILE\n" );

my $meta_container = $ENS_DBA->get_adaptor('MetaContainer');
my ($tax_id) = @{$meta_container->list_value_by_key('species.taxonomy_id')};
$tax_id || die( "\n[*DIE] No meta species.taxonomy_id in core DB!\n" );
my $genomedb_adaptor = $CMP_DBA->get_adaptor('GenomeDB');
my $mlss_adaptor     = $CMP_DBA->get_adaptor('MethodLinkSpeciesSet');
my $homology_adaptor = $CMP_DBA->get_adaptor('Homology');
my $member_adaptor   = $CMP_DBA->get_adaptor('Member');
my $tree_adaptor     = $CMP_DBA->get_adaptor('ProteinTree');

#---
#my $node = $tree_adaptor->fetch_node_by_node_id(699322);
#my $tree = $node->subroot;
#map{ warn('==> ', $_->node_id,' ',$_->parent->node_id) } 
#sort{$b->node_id<=>$a->node_id} @{$tree->get_all_nodes};
#die( $tree->node_id );
#---

my( $genome_db ) = ( grep{ $_->taxon_id == $tax_id } 
                     @{$genomedb_adaptor->fetch_all} );
$genome_db || die( "\n[*DIE] No taxon_id $tax_id in compara.\n" );

my $mlss = $mlss_adaptor->fetch_by_method_link_type_GenomeDBs
    ( $METHOD_LINK, [$genome_db] )
    || die( "\n[*DIE] No mlss for $METHOD_LINK and taxon $tax_id.\n" );

# Open the output file
open( OUT, "> $OUTFILE" ) 
    || die("\n[*DIE] Could not open $OUTFILE for write $!\n");
print OUT join( "\t", @FIELDS ). "\n";

#----------------------------------------------------------------------
# This is where the meat of the process is carried out
my %seen_genes;
my %bad_nodes; # Nodes that e.g. do not have leaves
my %homolog_with_no_tree;
my $cluster_count = 0;
MAIN:{
  
  my @homologs = @{$homology_adaptor
                       ->fetch_all_by_MethodLinkSpeciesSet($mlss)};
  #my @homologs = ($homology_adaptor->fetch_by_dbID(4));

  while( my $homology = shift @homologs  ){
    
    foreach my $gene_member( @{$homology->gene_list} ){
      next if $seen_genes{ $gene_member->stable_id }; # Process genes once
      my %record;

      my $tree = $tree_adaptor->fetch_by_Member_root_id($gene_member, 0);
      unless( $tree ){
        warn( "[WARN] No compara tree for ENSEMBLGENE "
              , $gene_member->stable_id );
        $homolog_with_no_tree{$gene_member->stable_id} ++;
        next;
      }

      my $node = $tree->get_leaf_by_Member($gene_member)
          || die( "[*DIE] No leaf for ENSEMBLGENE "
                  , $gene_member->stable_id );

      my $cluster_node = get_species_specific_cluster_node( $node );

      my @genes;
      foreach my $leaf( @{$cluster_node->get_all_leaves} ) {
        unless( $leaf->can('gene_member') ){ # Corrupt tree
          $bad_nodes{$leaf->node_id} ++;
          next;
        }
        my $gene_id = $leaf->gene_member->stable_id;
        $seen_genes{$gene_id} ++;
        push @genes, $gene_id;
      }
      next if scalar(@genes) < 2;
      
      $record{tree_node_id} = $tree->node_id;
      $record{tree_leaf_count} = $tree->num_leaves;
      $record{tree_taxon_name} = $tree->get_tagvalue('taxon_name');
      $record{cluster_node_id} = $cluster_node->node_id;
      $record{cluster_leaf_count} = scalar( @genes );
      $record{cluster_taxon_name} = $cluster_node->get_tagvalue('taxon_name');
      $record{cluster_genes} = join( ',', @genes );

      warn( "[INFO] Processed cluster: ", ++ $cluster_count );

      print OUT join
          ( "\t", 
            map{ $record{$_} } @FIELDS ) . "\n";
    
      #use Data::Dumper qw (Dumper);
      #warn Dumper( \%record );


    }


  }
} 

warn "[INFO] Total bad tree node ids: ". scalar(keys %bad_nodes) ."\n";
warn "[INFO] Bad node_ids: " . scalar(keys %bad_nodes) ."\n";
warn "[INFO] Homologs with no tree " . scalar (keys %homolog_with_no_tree);
warn "[INFO] Total genes: " . scalar(keys %seen_genes)."\n";

close OUT;

exit;

#======================================================================
# Takes a ProteinTree node, and recursively descends the tree until
# other species are seen, or the root node is reached. Returns the
# lowest node that contains all the same species
sub get_species_specific_cluster_node{
  my $node    = shift || return;
  my @leaves = @{$node->get_all_leaves};  
  my $species = $leaves[0]->genome_db->name;
  if( grep{ $_->can('genome_db') and ($_->genome_db->name ne $species) } 
      @leaves ){
    return ; # Mixed species
  }
  return get_species_specific_cluster_node( $node->parent ) || $node;
}

1;


