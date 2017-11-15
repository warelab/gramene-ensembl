#!/usr/bin/env perl

#paralog_clusters.pl dump_tree_id_by_species.out range > out

#e.g.
#range is the number of genes upstream to look for paralogs (e.g. 2 means that
#no more than 2 intervening non-paralogous genes would be allowed for a cluster.

use strict;

my ($tree_file, $range) = @ARGV;
$range++;

print_header();

my ($gene_fam, $gene_coord, $sorted_genes, $fam_taxon) = parse_tree_file($tree_file);

my $clust_cnt = 0;
my %clust_gene;  # e.g. $clust_gene{gid} = $clust_cnt;
my %tandem_cluster; #hash ref of array 

for my $chr (sort keys %$sorted_genes){
    my @gene_list = @{$sorted_genes->{$chr}};
    
    for my $i (0..$#gene_list){
        my ($gid, $s) = @{$gene_list[$i]};
	my $fam = $gene_fam->{$gid};
	
	my $positives = search_paralogs($i, $fam, \@gene_list); #array ref
	next unless $positives;
	
	for my $j (@$positives){
	    my ($pos_gid, $pos_s) = @{$gene_list[$j]};
	    my $pos_fam = $gene_fam->{$pos_gid};
            
	    if ($clust_gene{$gid}) {
	        my $current_clust = $clust_gene{$gid}; #cluster already defined
		next if $clust_gene{$pos_gid};
		$clust_gene{$pos_gid} = $current_clust;
		push @{$tandem_cluster{$current_clust}}, [$pos_gid, $j];
	    }
	    else {
	        $clust_cnt++;
		$clust_gene{$gid}     = $clust_cnt;
		$clust_gene{$pos_gid} = $clust_cnt;
		push @{$tandem_cluster{$clust_cnt}}, [$gid, $i]; 
		push @{$tandem_cluster{$clust_cnt}}, [$pos_gid,$j];
	    }
	}        
    }
}

for my $clust_id (sort {$a<=>$b} keys %tandem_cluster){
    my @members = @{$tandem_cluster{$clust_id}};
    
    for my $member (@members){
        my ($gid, $i) = @$member;
	my ($chr, $s, $e) = @{$gene_coord->{$gid}};
	my $root_id = $gene_fam->{$gid};
	my $taxon = $fam_taxon->{$root_id};
	print
	join(
	    "\t",
	    $clust_id,
            $root_id,
	    $gid,
	    $chr,
	    $s,
	    $e,
	    $i,
	    $taxon,
	), "\n";
    } 
}

sub search_paralogs {
    my ($i, $fam, $gene_list) = @_;
    my @positives;
    for my $j ($i + 1 .. $i + $range){
        last if $j >= $#{$gene_list}; #index of last element in array ref
	my ($gid, $s) = @{$gene_list->[$j]};
	push @positives, $j if $gene_fam->{$gid} eq $fam;
    }
    return \@positives;
}

sub parse_tree_file {
    my $file = shift;
    my %gene_fam;
    my ($unsorted, $sorted_genes);  #hash ref of arrays
    my %gene_coord;
    my %fam_taxon;

    open my $IN, "<$file" or die "can't open $file\n";
    while (<$IN>){
        next if /#/; #skip header
        chomp;
        my ($gid, $chr, $s, $e, $strand, $root_id, $tree_stable_id, $taxon) = split /\t/;
        next unless $root_id;
        #####
        #LOW CONFIDENCE GENE FILTERS
        #next if $taxon eq 'N'; # Actually Zea_mays (any of w22, b73, ph207)
        #next if $taxon eq 'Zea'; # Shared with parviglumis
	#next if $taxon eq 'Zea mays'; # used for any of the zea mays but means species-specific
	#next if $taxon eq 'Sorghum bicolor'; # sorghum specific
        ####
        $gene_fam{$gid} = $root_id;
	$fam_taxon{$root_id} = $taxon;
	push @{$unsorted->{$chr}}, [$gid, $s];
        $gene_coord{$gid} = [$chr, $s, $e];
    }
    
    #sort genes by coordinates
    for my $chr (sort keys %$unsorted){
        push @{$sorted_genes->{$chr}}, sort {$a->[1]<=>$b->[1]} @{$unsorted->{$chr}};
    }
    
    return (\%gene_fam,\%gene_coord, $sorted_genes, \%fam_taxon);
}

sub print_header {
    print
    join
    ("\t",
     'clust_id',
     'root_id',
     'gid',
     'chr',
     'start',
     'end',
     'index',
     'root_taxon',
    ), "\n";
}
