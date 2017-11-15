#!/usr/bin/env perl

#add_coords_to_syn_table ort_file syn_table > out

use strict;


=stub

my %sp = ('4577' => 'Zm',
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

my ($ort_file, $align_file, $species_meta, $taxid2sp_file) = @ARGV;

#my $taxid2sp_file = shift;


my %sp; #
if( -e $taxid2sp_file){
	%sp = %{get_taxid2sp_hash( $taxid2sp_file )};
}



if ($species_meta) {
	foreach my $meta (split /,/, $species_meta) {
		$meta=~/(\w+)=(\w+)/;
		if ($1 && $2) {
			$sp{$1}=$2;
		}
	}
}

my ($coord) = get_coords($ort_file);

my ($block1, $block2, $ord) = get_dag_info($align_file);

print_header();

for my $block_id (sort keys %$block1){
    my (undef, undef, $strand) = split "_", $block_id;
    $strand =~ s/\.\d.*//;
    my $str = $strand eq 'F' ? '+' : '-';
	
    my @gid1 = @{$block1->{$block_id}};
    my @gid2 = @{$block2->{$block_id}};
    
    my ($ref1, $s1, $e1, $ord1_s, $ord1_e, $tax1) = block_coords(\@gid1);
    my ($ref2, $s2, $e2, $ord2_s, $ord2_e, $tax2) = block_coords(\@gid2);
    
    my $range1 = '(' . $ord1_s . '..' . $ord1_e . ')';
    my $range2 = $str eq '+' ? '(' . $ord2_s . '..' . $ord2_e . ')'
                             : '(' . $ord2_s . '..' . $ord2_e . ')';
    
    print
    join("\t",
        ($taxid2sp_file ? $sp{$tax2} : $tax2) . $ref2,
        $ref2,
	$s2,
	$e2,
	'+',
#	$ord1_s,
#	$ord1_e,
	$range2,
	$tax2,
	($taxid2sp_file ? $sp{$tax1} : $tax1)  . $ref1,
	$ref1,
	$s1,
	$e1,
	$str,
#	$ord2_s,
#	$ord2_e,
	$range1,
	$tax1,
	$block_id,
    ), "\n";
    
}

sub block_coords {
    my $gids = shift;
    my $block_ref;
    my $block_tax;
    my @coords;
    my @ords;
    
    for my $gid (@$gids){
        my ($ref, $s, $e, $tax) = @{$coord->{$gid}};
	push @coords, $s, $e;
	
	push @ords, $ord->{$gid};
	
	$block_ref = $ref;
	$block_tax = $tax;
    }
    
    my @sorted = sort {$a<=>$b} @coords;
    my $s = $sorted[0];
    my $e = $sorted[-1];
    
    my @sorted_ord = sort {$a<=>$b} @ords;
    my $ord_s = $sorted_ord[0];
    my $ord_e = $sorted_ord[-1];
    
    return ($block_ref, $s, $e, $ord_s, $ord_e, $block_tax);
    
}


sub get_dag_info {
    my $align_file = shift;
    my (%block1, %block2, %ord);
    my $block_id;
    
    open my $IN, "<$align_file" or die "can't open $align_file\n";
    while(<$IN>){
        chomp;
        my @f=split;
        if ($f[0] eq '##'){  #block header
            my ($aln_num, $score, $pair_cnt, $refA, $refB) 
                = ($f[-8], $f[-5], $f[-1], $f[2], $f[4]);
            $aln_num =~ s/\#//;
            $pair_cnt =~ s/\)://;
            my $strand = $f[5] eq '(reverse)' ? 'R' : 'F';
            $block_id = $refB . '_' . $refA . '_' . 
                        $strand . '.' . $aln_num . 
                        '(' . $pair_cnt . ')';  #eg: Zm1_Sb1_F.13(7)
        }
        else {
	     my ($gid1, $gid2) = ($f[5], $f[1]);
	     my ($ord1, $ord2) = ($f[6], $f[2]);
	     
	     push @{$block1{$block_id}}, $gid1;
	     push @{$block2{$block_id}}, $gid2;
	     
	     $ord{$gid1} = $ord1;
	     $ord{$gid2} = $ord2;
	     
#            my ($ref1, $ref2) = ($_[4], $_[0]);
            
#            push @{$block{$gid1}->{$gid2}}, $block_id;
#            $aln1_2{$gid1}->{$gid2} = 1;
#            push @{$dagorder1{$ref1}}, [$gid1, $ord1,] unless $daggene1{$gid1};
#            $daggene1{$gid1} = [$ref1, $ord1,];
        }
    }
    
    return (\%block1, \%block2, \%ord);

}


sub get_coords {
    my $file = shift;
    my %coord;
    
    open my $IN, "<$file" or die "can't open $file\n";
    while(<$IN>){
       chomp;
       next if /gene_stable/; #skip header
       my (
           $gid1,
           $ref1,
           $s1,
           $e1,
           $tax1,
           $gid2,
           $ref2,
           $s2,
           $e2,
           $relationship,
           $tax2,
           ) = split /\t/;
       
       $coord{$gid1} = [$ref1, $s1, $e1, $tax1];
       $coord{$gid2} = [$ref2, $s2, $e2, $tax2];
    }
    
    return (\%coord);
}

sub print_header {
    print
    join("\t",
        'sp1_chr',
	'sp1_ref',
	'sp1_s',
	'sp1_e',
	'sp1_str',
	'sp1_range',
	'sp1_taxID',
	'sp2_chr',
	'sp2_ref',
	'sp2_s',
	'sp2_e',
	'sp2_str',
	'sp2_range',
	'sp2_taxID',
	'blockID',	
    ), "\n";
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
                        print STDERR "Weix:Unrecoganized species $_\n";
                }
        }

        return \%tax2sp_hash;
}

