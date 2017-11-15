#!/usr/bin/perl

#combine_tables.pl
#Joshua Stein
#Modified 20061205

use strict;
use Getopt::Long;

if (@ARGV < 4) {die "\nusage: combine_tables.pl file1 index1 file2 index2 [-h] > out
\nColumns in file1 are combined with file2 based on common
key as specified by column numbers in index1 and index2.\n
Option -h combines headers without requiring a common key.
Do not use this option unless both files have column headers.\n\n" };

my $header;
GetOptions('header' => \$header);
my($file1, $index1, $file2, $index2) = @ARGV;


if ($header){
    process_w_header();
} else {
    process_wo_header();
}    

sub process_w_header{
    my %hash;
    #my %hash2;
    my $max_field_count;
    $index1 = $index1 - 1;
    $index2 = $index2 - 1;
    my $line_count1;
    my $line_count2;
    my $header1;
    my $header2;
    my @lines;
    
    open(TAB, "<$file2") || die "Cannot open $file2\n";
    while (<TAB>) {
        chomp;
        $line_count2++;
	if ($line_count2 == 1){
	    $header2 = $_;
	    #print $_, "\t";
	    next;
	}
        my @fields = split (/\t/);
        my $field_count = @fields;
        $max_field_count = $field_count unless $field_count < $max_field_count;
        my $key = $fields[$index2];
        push (@{$hash{$key}}, $_ );
    }
    close (TAB);

    my $filler = "\t" x $max_field_count;

    open(TABTWO, "<$file1") || die "Cannot open $file1\n";
    while (<TABTWO>) {
        chomp;
	$line_count1++;
	if ($line_count1 == 1){
	    $header1 = $_;
	    #print $_, "\n";
	    next;
	}
	
        push @lines, $_;
        #split (/\t/);
        #my $key = $_[$index1];
        #push (@{$hash1{$key}}, $_ );
    }
    close (TABTWO);
    
    print "$header1\t$header2\n";
    
    for my $line (@lines){
        split(/\t/, $line);
	my $key = $_[$index1];    
        if ( exists( $hash{$key} ) ) {
            foreach my $value ( @{ $hash{$key} } ) {
	        print $line,"\t", $value, "\n";
            }
        } else {
            print $line, $filler,"\n";
	}
    }#end for
} #end process_wo_header

sub process_wo_header{
    #This is the original script:
    
    my %hash;
    my $max_field_count;
    $index1 = $index1 - 1;
    $index2 = $index2 - 1;
    open(TAB, "<$file2") || die "Cannot open $file2\n";
    while (<TAB>) {
        chomp;
        my @fields = split (/\t/);
        my $field_count = @fields;
        $max_field_count = $field_count unless $field_count < $max_field_count;
        my $key = $fields[$index2];
        push (@{ $hash{$key} }, $_ );
    }
    close (TAB);

    my $filler = "\tNULL" x $max_field_count;

    open(TABTWO, "<$file1") || die "Cannot open $file1\n";
    while (<TABTWO>) {
        chomp;
        split (/\t/);
        my $key = $_[$index1];
        if ( exists( $hash{$key} ) ) {
            foreach my $value ( @{ $hash{$key} } ) {
	        print "$_\t$value\n";
            }
        } else {
            print "$_$filler\n";
	}
    }
    close (TABTWO);
} #end process_wo_header
