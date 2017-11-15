#!/usr/bin/env perl

#ort2matchlist file.ort > file.matchlist

use strict;

my ($ort_file) = @ARGV;

parse_ort_file( $ort_file );

sub parse_ort_file {
    my $file = shift;
    my (%gene_data, %ort, %gene_list, %tax1, %tax2);
    
    open IN, "<$file" or die "cannot open file $file:$!\n";
    while (<IN>){
        chomp;
        next if /gene_stable_id|gid1/; #skip header
        my ($gid1,$gid2,$coord1,$coord2,$sp1,$sp2,$tax1,$tax2,$hom_type,$tree_root_id)=split/\t/;
	my ($chr1,$s1,$e1)=split(/:/,$coord1);
	my ($chr2,$s2,$e2)=split(/:/,$coord2);

	my $new_ort_line = join "\t", ($gid1,$chr1,$s1,$e1,$tax1,$gid2,$chr2,$s2,$e2,$hom_type,$tax2);
	print "$new_ort_line\n";
    }
    close IN;
    
}
