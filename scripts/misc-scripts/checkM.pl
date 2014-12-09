#!/bin/env perl

use Bio::SeqIO;

my $f = shift;

my $seqI = Bio::SeqIO->new(-file => $f, -format => 'fasta');

my $count;

while (my $s=$seqI->next_seq){
 	if ($s->seq =~ /^M/){
		next;
	}else{
		my $ID=$s->display_id;
		print "$ID not start with M\n";
		$count++;
	}	
}

print "In total $count transcripts do not start with M\n";

