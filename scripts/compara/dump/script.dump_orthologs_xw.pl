#! /usr/bin/env perl

use strict;
use warnings;
use Data::Dumper qw(Dumper);
use File::Basename;
use FindBin qw( $Bin );
use Pod::Usage;
use Getopt::Long;
use IO::File;
use Bio::EnsEMBL::Registry;

#our $ENS_DBA;
#our $CMP_DBA;

my( $sp1, $sp2, $reg,);
GetOptions ( 
    "species=s"          => \$sp1,
    "ensembl_registry=s" => \$reg,
    "other=s"            => \$sp2,
    );

print_header();
# Load the ensembl file and adapters
Bio::EnsEMBL::Registry->load_all( $reg );

#adaptor for homologies
my $gene_member_adaptor = Bio::EnsEMBL::Registry->get_adaptor('compara', 'compara', 'GeneMember');
my $homology_adaptor = Bio::EnsEMBL::Registry->get_adaptor('compara', 'compara', 'Homology');

# Get all genes from refernce species
my $gene_adaptor = Bio::EnsEMBL::Registry->get_adaptor( $sp1, 'core', 'Gene' );
my $gene_id_arrayref = $gene_adaptor->list_stable_ids;

while( my $gid = shift @$gene_id_arrayref ){ #a gene stable id
    #print "$gid\n";
    next unless $gid;
    my $gm1 = $gene_member_adaptor->fetch_by_stable_id($gid); #gene_member
    next unless $gm1;
    my $homologies = $homology_adaptor->fetch_all_by_Member($gm1, -TARGET_SPECIES => $sp2); #this should filter for target but doesn't work
    next unless $homologies;

    for my $homology (@{$homologies}) {
	my $hom_type = $homology->description;
	next unless $hom_type =~ /ortholog/;

	for my $gm2 ( @{$homology->gene_list()} ){ #gene_member
	    next unless $gm2->genome_db()->name eq $sp2;

	    my $coord1 = join(":", $gm1->dnafrag()->name(), $gm1->dnafrag_start, $gm1->dnafrag_end, $gm1->dnafrag_strand,);
	    my $coord2 = join(":", $gm2->dnafrag()->name(), $gm2->dnafrag_start, $gm2->dnafrag_end, $gm2->dnafrag_strand,);

	    print join("\t", 
		       $gid, 
		       $gm2->stable_id,
		       $coord1,
		       $coord2,
		       $sp1,
		       $sp2,
		       $gm1->taxon_id . ucfirst(join '', map { (substr $_, 0, 1).(substr $_, -1)} (split '_', $sp1)), #ucfirst (join '', map{ substr($_, 0, 1) } (split '_', $sp1)),
		       $gm2->taxon_id . ucfirst (join '', map{ substr($_, 0, 1).(substr $_, -1) } (split '_', $sp2)),
		       $hom_type,
		       $homology->gene_tree()->root_id(),
		), "\n";
	}
    }
}

sub print_header {
    my @header = qw/gid1 gid2 coord1 coord2 sp1 sp2 taxid1 taxid2 hom_type tree_root_id/;
    print join("\t", @header), "\n";
}
