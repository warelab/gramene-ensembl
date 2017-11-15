#! /usr/bin/env perl

#./script.correct_split.pl result.count_splits.out tree_id.zea_mays > tree_id.corrected.zea_mays

use strict;

my ($split_file, $tree_id_file) = @ARGV;

my $splits = get_splits($split_file);
my %seen;
open my $IN, "<$tree_id_file" or die;
while(<$IN>){
    print $_ and next if /#/;
    chomp;
    my ($gid, $ref_name, $start, $end, $strand, 
	$root_node_id, $tree_stable_id, $root_taxon) = split /\t/;
    if ($splits->{$gid}){
	my ($fused_gid, $fused_chr, $fused_s, $fused_e) = @{$splits->{$gid}};
	next if $seen{$fused_gid};
	$seen{$fused_gid} = 1;
	print join("\t", $fused_gid, $ref_name, $fused_s, $fused_e, $strand, $root_node_id, $tree_stable_id, $root_taxon),"\n";
    }
    else {
	print $_, "\n";
    }
}

sub get_splits {
    my $file = shift;
    my %splits;
    open my $IN, "<$file" or die;
    while(<$IN>){
	next if /node_taxon/; #skip header
	chomp;
	my ($split_gid, $split_chr, $split_s, $split_e, $split_strand, 
	    $fused_gid, $fused_chr, $fused_s, $fused_e, $fused_strand,
	    $species, $tree_id, $node_id) = split /\t/;
	$splits{$split_gid} = [$fused_gid, $fused_chr, $fused_s, $fused_e];
    }
    return \%splits;
}
