#!/usr/bin/env perl

#ort2matchlist file.ort > file.matchlist

use strict;

my $ort_file = shift;
my $species_meta = shift;
my $taxid2sp_file = shift;

#warn ("weix debug get arg ort_file=$ort_file, species_meta=$species_meta, taxid2sp_file=$taxid2sp_file\n");
my %sp;
%sp = %{get_taxid2sp_hash( $taxid2sp_file )} if $taxid2sp_file;

for my $k (keys %sp){
	print STDERR "$k => $sp{$k}\n";
} 

=stub

('4577' => 'Zm',
          '4558' => 'Sb',
          '39947' => 'Os',
	  '15368' => 'Bd',
	  '3702'	=> 'At',
	  '59689'	=> 'Al',
	  '29760'	=> 'Vv',
	  '3694'	=> 'Pt',
	  '4538'	=> 'Og',
	  '4533'	=> 'Ob',
	  '4555'	=> 'Si',
	  '39946'	=> 'OsIG',
	  '112509'	=> 'Hv',
	  '214687'	=> 'Ma',
	  '4572'	=> 'Tu',
	  '37682'	=> 'At',
	  '65489' => 'Ob',
	  '4533' => 'Ob',
	  '4538' => 'Og',
	  '40148' => 'Og',
	  '39946' => 'Oi',
	  '40149' => 'Om',
	  '4536' => 'On',
	  '4537' => 'Op',
	  '4529' => 'Or',
	  '39947' => 'OsJG',
	  '77586' => 'Lp',
          );
=cut


if ($species_meta) {
	foreach my $meta (split /,/, $species_meta) {
		$meta=~/(\w+)=(\w+)/;
		if ($1 && $2) {
			$sp{$1}=$2;
		}
	}
}

#(my $sp1_sp2 = $ort_file) =~ s/\.ort//;
#$sp1_sp2 =~ s/.*\///;  #strip off leading directory
#my ($sp1, $sp2) = split '_', $sp1_sp2;


my ($gene_data, $ort, $tax1, $tax2) = parse_ort_file($ort_file);

for my $gid1 (keys %{$ort}){
    my $tax1 = $tax1->{$gid1};
    my ($ref1, $idx1) = @{$gene_data->{$tax1}->{$gid1}};
    for my $gid2 (keys %{$ort->{$gid1}}){
        my $tax2 = $tax2->{$gid2}; 
	my ($ref2, $idx2) = @{$gene_data->{$tax2}->{$gid2}};
	if ( $ref1 gt $ref2) {
	print
	join("\t",
	    $ref2,
        $gid2,
	    $idx2,
	    $idx2,
	    $ref1,
	    $gid1,
	    $idx1,
	    $idx1,
	    0,
	), "\n";
	} else {
	print
	join("\t",
	    $ref1,
        $gid1,
	    $idx1,
	    $idx1,
	    $ref2,
	    $gid2,
	    $idx2,
	    $idx2,
	    0,
	), "\n";
	}
    }
}

sub parse_ort_file {
    my $file = shift;
    my (%gene_data, %ort, %gene_list, %tax1, %tax2);
    
    open IN, "<$file" or die "cannot open file $file:$!\n";
    while (<IN>){
        chomp;
        next if /gene_stable_id/; #skip header
        my ($gid1,    #gene id of first species
            $ref1,    #Ref chr or BAC accession if sp1 is maize (4577)
            $s1,      #Start of gid1
            $e1,      #End of gid1
            $tax1,     #species 1 id (e.g. 4577Zm maize or 39947Os rice japonica)
            $gid2,    #Same info for species 2
            $ref2,
            $s2,
            $e2,
            $hom_relationship,  #e.g. ortholog_one2many
            $tax2, 
            ) = split /\t/;
	    
        my $sp1 = $taxid2sp_file ? $sp{$tax1} : $tax1;
	my $sp2 = $taxid2sp_file ? $sp{$tax2} : $tax2;
	
	$ref1 = $sp1 . $ref1;
	$ref2 = $sp2 . $ref2;
	
        push @{$gene_list{$tax1}->{$ref1}}, [$gid1, $s1, $e1,]
	    unless $tax1{$gid1};  #push this gene only once
	                          #otherwise will get gap in order
	push @{$gene_list{$tax2}->{$ref2}}, [$gid2, $s2, $e2,]
	    unless $tax2{$gid2};  #push this gene only once
	                          #otherwise will get gap in order;
        $ort{$gid1}->{$gid2} = 1;
	$tax1{$gid1} = $tax1;
	$tax2{$gid2} = $tax2;
    }
    close IN;
    
    for my $tax (keys %gene_list){
        for my $ref (keys %{$gene_list{$tax}} ) {
            my @sorted = sort {$a->[1]<=>$b->[1]} @{$gene_list{$tax}->{$ref}};	    
	    for my $idx (0..$#sorted){
	        my $gid = $sorted[$idx]->[0];
		my $s   = $sorted[$idx]->[1];
		my $e   = $sorted[$idx]->[2];
	        $gene_data{$tax}->{$gid} = [$ref, $idx, $s, $e]; 
#		print "$sp\t$chr\t$idx\t$gid\n";
	    } 
        }
    } 
    
    return (\%gene_data, \%ort, \%tax1, \%tax2,);
}


sub get_taxid2sp_hash{

	my $file = shift;
	my %tax2sp_hash;

	open my $fh, "<$file" or die "cannot open file $file:$!\n";
    	while (<$fh>){
        	chomp;
		if(/^(\d+)\s+(.+)$/){
			$tax2sp_hash{$1} = ucfirst (join '', map{ substr($_, 0, 1) } (split '_', $2)); 	
		}else{
			print STDERR "Weix: Unrecoganized species $_\n";
		}
	}

	return \%tax2sp_hash;
}
