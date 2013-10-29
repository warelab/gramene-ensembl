#!/usr/bin/env perl

package main;
our $SEE;

package Nuc_translator;

use strict;
require Exporter;
our @ISA = qw (Exporter);
our @EXPORT = qw (translate_sequence get_protein reverse_complement);

our $currentCode;
our %codon_table;


## See http://golgi.harvard.edu/biolinks/gencode.html
my %SUPPORTED_GENETIC_CODES = ( universal => 1,
				Euplotes => 1,
				Tetrahymena => 1,
				Candida => 1,
				Acetabularia => 1);


sub translate_sequence {
    &_init_codon_table() unless $currentCode;
    my ($sequence, $frame) = @_;

    $sequence = uc ($sequence);
    $sequence =~ tr/T/U/;
    my $seq_length = length ($sequence);
    unless ($frame > 0 and $frame < 4) { $frame = 1;}
    my $start_point = $frame - 1;
    my $protein_sequence;
    for (my $i = $start_point; $i < $seq_length; $i+=3) {
	my $codon = substr($sequence, $i, 3);
	my $amino_acid;
	if (exists($codon_table{$codon})) {
	    $amino_acid = $codon_table{$codon};
	} else {
	    if (length($codon) == 3) {
		$amino_acid = 'X';
	    } else {
		$amino_acid = "";
	    }
	}
	$protein_sequence .= $amino_acid;
    }
    return($protein_sequence);
}

sub get_protein {
    my ($sequence) = @_;
    
    ## Assume frame 1 unless multiple stops appear.
    my $least_stops = undef();
    my $least_stop_prot_seq = "";
    foreach my $forward_frame (1, 2, 3) {
	my $protein = &translate_sequence($sequence, $forward_frame);
	my $num_stops = &count_stops_in_prot_seq($protein);
	if ($num_stops == 0) {
	    return ($protein);
	} else {
	    if (!defined($least_stops)) {
		#initialize data
		$least_stops = $num_stops;
		$least_stop_prot_seq = $protein;
	    } elsif ($num_stops < $least_stops) {
		$least_stops = $num_stops;
		$least_stop_prot_seq = $protein;
	    } else {
		#keeping original $num_stops and $least_stop_prot_seq
	    }
	}
    }
    return ($least_stop_prot_seq);
}

sub reverse_complement {
    my($s) = @_;
    my ($rc);
    $rc = reverse ($s);
    $rc =~tr/ACGTacgtyrkmYRKM/TGCAtgcarymkRYMK/;
    return($rc);
}


####
sub count_stops_in_prot_seq {
    my ($prot_seq) = @_;
    chop $prot_seq; #remove trailing stop.
    my $stop_num = 0;
    while ($prot_seq =~ /\*/g) {
	$stop_num++;
    } 
    return ($stop_num);
}



####
sub use_specified_genetic_code {

    my ($special_code) = @_;
    print STDERR "using special genetic code $special_code\n" if $SEE;
    unless ($SUPPORTED_GENETIC_CODES{$special_code}) {
	die "Sorry, $special_code is not currently supported or recognized.\n";
    }
    _init_codon_table(); ## Restore default universal code.  Others are variations on this.
    $currentCode = $special_code;
    
    if ($special_code eq "Euplotes") {
	$codon_table{UGA} = "C";
    } 

    elsif ($special_code eq "Tetrahymena" || $special_code eq "Acetabularia") {
	$codon_table{UAA} = "Q";
	$codon_table{UAG} = "Q";
    }
    
    elsif ($special_code eq "Candida") {
	$codon_table{CUG} = "S";
    }

}


####
sub get_stop_codons {
    my @stop_codons;
    foreach my $codon (keys %codon_table) {
	if ($codon_table{$codon} eq '*') {
	    push (@stop_codons, $codon);
	}
    }
    foreach my $codon (@stop_codons) {
	$codon =~ tr/U/T/;
    }
    return (@stop_codons);
}


sub _init_codon_table {
    print STDERR "initing codon table.\n" if $SEE;
    ## Set to Universal Genetic Code
    $currentCode = "universal";

    %codon_table = (    UUU => 'F',
			UUC => 'F',
			UUA => 'L',
			UUG => 'L',
			
			CUU => 'L',
			CUC => 'L',
			CUA => 'L',
			CUG => 'L',
			
			AUU => 'I',
			AUC => 'I',
			AUA => 'I',
			AUG => 'M',
			
			GUU => 'V',
			GUC => 'V',
			GUA => 'V',
			GUG => 'V',
			
			UCU => 'S',
			UCC => 'S',
			UCA => 'S',
			UCG => 'S',
			
			CCU => 'P',
			CCC => 'P',
			CCA => 'P',
			CCG => 'P',
			
			ACU => 'T',
			ACC => 'T',
			ACA => 'T',
			ACG => 'T',
			
			GCU => 'A',
			GCC => 'A',
			GCA => 'A',
			GCG => 'A',
			
			UAU => 'Y',
			UAC => 'Y',
			UAA => '*',
			UAG => '*',
			
			CAU => 'H',
			CAC => 'H',
			CAA => 'Q',
			CAG => 'Q',
			
			AAU => 'N',
			AAC => 'N',
			AAA => 'K',
			AAG => 'K',
			
			GAU => 'D',
			GAC => 'D',
			GAA => 'E',
			GAG => 'E',
			
			UGU => 'C',
			UGC => 'C',
			UGA => '*',
			UGG => 'W',
			
			CGU => 'R',
			CGC => 'R',
			CGA => 'R',
			CGG => 'R',
			
			AGU => 'S',
			AGC => 'S',
			AGA => 'R',
			AGG => 'R',
			
			GGU => 'G',
			GGC => 'G',
			GGA => 'G',
			GGG => 'G'    
			
			);
}



1; #end of module
