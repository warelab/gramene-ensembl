#!/usr/bin/env perl
#
##get_g2pcid.pl pc.3.species > out
#
##e.g.
#parse the paralog cluster output file and get mapping of gene to paralog cluster ID

use strict;

my ($pc_file) = @ARGV;
#
#
#cat pc.3.w22.out | grep -v 'clust_id' | cut -f1,3,8 | perl -e 'my %clust; while(<>){chomp; my($cl,$gid)=split/\t/; push @{$clust{$cl}}, $gid} for my $cl (keys %clust){my @gids = @{$clust{$cl}}; my $cnt = scalar @gids; my $id = "pc$cl.w22.$cnt"; for my $gid (@gids){print join("\t", $gid,$id),"\n";} }' > gid_pcID.w22

my %clust;
my $fh;

my $sp;
if( $pc_file =~ /.*\.(\w+)/){
	$sp = $1;
}
if( $sp =~ /([a-z0-9])[a-z0-9]*_(\w+)/){
	$sp = uc($1).$2;
	$sp =~ s/_//g;
}

open $fh, $pc_file or die "Cannot open $pc_file to read";

while ( <$fh> ){
	next if /clust_id/i;
	chomp;

	my @fields = split /\t/;
	my $cluster_id = $fields[0];
	my $gid        = $fields[2];
	
	push @{$clust{$cluster_id}}, $gid;
}

close $fh;

for my $cl (keys %clust){

	my @gids = @{$clust{$cl}}; 
	my $cnt = scalar @gids; 
	my $id = "pc$cl.$sp.$cnt"; 
	
	for my $gid (@gids){
		print join("\t", $gid,$id),"\n";
	} 

}
