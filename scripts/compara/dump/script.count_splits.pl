#! /usr/bin/env perl

#usage: count_splits.pl -e ensembl.registry > out
#script loops over trees, then loops over nodes, then filter for 'gene_split', then collect leaves, then get gene info
#'node_id' associates split genes with one another

use strict;
#use warnings;
use Getopt::Long;
use Bio::EnsEMBL::Registry;

my( $reg,);
GetOptions ( 
    "ensembl_registry=s" => \$reg,
    );

# Load the ensembl file and adapters
Bio::EnsEMBL::Registry->load_all( $reg );

#need a tree adaptor
my $tree_adaptor = Bio::EnsEMBL::Registry->get_adaptor('compara', 'compara', 'GeneTree');

#collect protein-coding trees
my @trees = @{$tree_adaptor->fetch_all(-tree_type     => 'tree',
				       -member_type   => 'protein',
				       -clusterset_id => 'default',)};

my @header = (qw/species gid coord node_type tree_id node_taxon node_id/);
print join("\t", @header), "\n";

#loop over trees, then loop over nodes, then filter for 'gene_split', then collect leaves, then get gene info 
while (my $tree = shift @trees){
    my $tree_taxon = $tree->root->taxonomy_level();
    my $tree_id    = $tree->stable_id();
    my $node_arrayref = $tree->get_all_nodes; #doesn't work: my $node_arrayref = $tree->get_all_nodes_by_tag_value('node_type' => 'gene_split');
    collect_split_node_data($node_arrayref, $tree_id, $tree_taxon);
    $tree->release_tree();
}

sub collect_split_node_data {
    my ($nodes, $tree_id, $tree_taxon) = @_; #node_arrayref;  tree_id is the root node id of the tree
    my %seen;
    for my $node( @$nodes ){ #Bio::EnsEMBL::Compara::GeneTreeNode object 
	next if $node->is_leaf;
	my $node_type = $node->node_type(); 
	next unless $node_type eq 'gene_split';
        my $node_taxon_name = $node->taxonomy_level(); # e.g. Magnoliophyta
	my $node_id = $node->node_id; #id of gene_split node
	my @leaves = @{$node->get_all_leaves() }; #arrayref Bio::EnsEMBL::Compara::NestedSet

	for my $leaf (@leaves){ #leaf is a Bio::EnsEMBL::Compara::SeqMember (stable_id is a protein or transcript)
	    my $gene = $leaf->gene_member; #converts to a Bio::EnsEMBL::Compara::SeqMember::gene_member
	    my $gid = $gene->stable_id();
	    next if $seen{$gid};  #if gene_split node is not terminal to the leaves then the same gene will get revisited in child nodes which we don't want 
	    $seen{$gid} = 1;
	    my $coord = join(":", $gene->dnafrag()->name(), $gene->dnafrag_start, $gene->dnafrag_end, $gene->dnafrag_strand,);
	    print join("\t", 
		       $leaf->genome_db()->name, #species e.g. zea_mays
		       $gid,
                       $coord,
		       $node_type, #should always be 'gene_split'
		       $tree_id,   #stable id, e.g. EPlGT00820000103153
		       $node_taxon_name,
		       $node_id, #gene models with same node_id are split from the same gene
		      # $tree_taxon,
		), "\n";
	}
    }
}
