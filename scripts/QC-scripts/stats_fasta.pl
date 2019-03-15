#! /usr/bin/env perl

use strict;
use Bio::SeqIO;
use Statistics::Descriptive;

if (@ARGV < 1) { die "Usage: seq_stats.pl input_1.fa [input_2.fa] ...[input_n.fa]\n"; }

#print header
print join(
    "\t",
    #'file',
    'count',
    'sum_len',
    'N50',
    'min_len',
    'max_len',
    'med_len',
    'ave_len',
    'sd_len',
    'file',
), "\n";

#Loop over each file, print out results
for my $input ( @ARGV ) {

    my $lengths = get_lengths($input); #arrayref of sequences lengths
    
    my ($cnt, $sum, $min, $max, $med, $ave, $sd) = get_stats($lengths);
    
    my $n50 = get_n50($lengths, $sum);
    
    print join(
        "\t",
        #$input,
        $cnt,
        $sum,
	$n50,
        $min,
        $max,
	$med,
	$ave,
	$sd,
	$input
    ), "\n";
} #end for


sub get_lengths{
    my $input = shift;
    my @lengths;
    my $in  = Bio::SeqIO-> new ( -file   => $input,
                                 -format => 'fasta' );
			     
    while (my $record = $in->next_seq()) {
         my $id = $record->id();
         my $seq = $record->seq();
         my $seq_len = length $seq;
         push @lengths, $seq_len;
    }
    return \@lengths;
}


sub get_stats{
    my $lengths = shift;
    
    my $stat = Statistics::Descriptive::Full->new();
    $stat->add_data( @{$lengths} );
    
    my $cnt = $stat->count();    
    my $sum = $stat->sum();
    my $min = $stat->min();
    my $max = $stat->max();
    my $med = $stat->median();
    
    my $ave = $stat->mean();
    $ave = sprintf("%.0f", $ave);
    
    my $sd   = $stat->standard_deviation();
    $sd = sprintf("%.0f", $sd);
    
    return ($cnt, $sum, $min, $max, $med, $ave, $sd);
}

sub get_n50{
    my ($lengths, $sum) = @_;
    my @lengths = sort {$a<=>$b} @{$lengths};
    my $subtotal;
    #my $n50 = 0;
    while (my $length = shift @lengths){
        $subtotal += $length;
	if ($subtotal >= $sum/2){
	    return $length;
	}
	#print "$length\t$subtotal\t$n50\n";	
    }
    
    #print join("\n", @lengths), "\n";
}
