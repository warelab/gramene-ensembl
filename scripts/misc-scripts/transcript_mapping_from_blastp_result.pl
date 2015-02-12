#!/bin/env perl

#parse blastp result and identify the best matching v1 transcrpt
#to v3 transcript, and  v1 gene for each v3 gene


use strict;
use Getopt::Long;
use Data::Dumper;
use Pod::Usage;
use File::Basename;
use Cwd 'abs_path';
use vars qw[ $VERSION ];
use File::Path qw(make_path remove_tree);
use Number::Interval;


$VERSION = sprintf "%d.%02d", q$Revision: 1.10 $ =~ /(\d+)\.(\d+)/;
pod2usage(0) unless $ARGV[0];

my ($_help, $_version, $_debug, $infile, $outfile);
GetOptions(
	   'h|help'        => \$_help,
	   'v|version'     => \$_version,
	   'infile=s'      => \$infile,
           'outfile=s'     => \$outfile,
           'debug'         => \$_debug,
	   ) or die;

pod2usage(2) if $_help;
if ( $_version ) {
    my $prog = basename( $0 );
    print "$prog v$VERSION\n";
    exit(0);
}

our $keys = {
    1 => 'qid',
    2 => 'sid',
    3 => 'qlen',
    4 => 'slen',
    5 => 'evalue',
    6 => 'score',
    7 => 'qstart',
    8 => 'qend',
    9 => 'sstart',
    10 => 'send',
    11 => 'gapopen',
    12 => 'gaps',
    13 => 'pident',
    14 => 'length',
}; 

# 1. query_trpt_stable_id|chr:chr_start:chr_end:chr_strand:translation_stable_id:gene_stable_id   OMERI12G00730.1|12:577159:589003:1:OMERI12G00730.1:OMERI12G00730
# 2. subject_trpt_stable_id|chr:chr_start:chr_end:chr_strand:translation_stable_id:gene_stable_id OMERI12G01440.1|12:1057633:1068367:-1|OMERI12G01440.1|OMERI12G01440   # 3. qlen  96
# 4. slen      605
# 5. evalue     0.66
# 6. score    63
# 7. qstart      7
# 8. qend       34
# 9. sstart      571
# 10. send     598
# 11. gapopen     0
# 12. gaps       0
# 13. pident       39.29
# 14. length   28


my $ifh;
if ($infile =~ /\.gz$/){
    open $ifh, "gunzip -c $infile |" or die "Cannot open $infile to read";
}else{
    open $ifh, $infile or die "Cannot open $infile to read";
}

my %qid_hits_hash;
while (my $ahit_line = <$ifh> ){
    print "$ahit_line" if $_debug;
    chomp $ahit_line;
    my @fields = split ' ', $ahit_line;
    
    my $qID = shift @fields;
    my $i=2;
    my %hash =  map { $keys->{$i++} => $_ } @fields;
    map {print "$_, $hash{$_}\n" } keys %hash if $_debug;
    push @{$qid_hits_hash{$qID}}, \%hash;
}

close $ifh;



our %qid2tidmapping; #populated by funciton id_mapping
our %qidnomapping; #populated by funciton id_mapping

foreach my $trptID ( sort keys %qid_hits_hash){
    id_mapping($trptID, $qid_hits_hash{$trptID});
    print "processing $trptID\n" if $_debug;
    
}

#exit;

print_id_mapping();

print_id_nomapping();

#====================================
# Function to pick the best match
# target (usually older version) gene
# for query (new eversion) gene

sub id_mapping{
    
    my $trpt_id = shift;
    my $hits_arrayref = shift;

    my @sorted_hits = sort {$a->{evalue} <=> $b->{evalue}, 
			    $b->{pident} <=> $a->{pident},
			    $a->{gapopen} <=> $b->{gapopen},
                           } @{$hits_arrayref};

 
    my $overlap = 0;
    my $best_hit = undef;
    foreach my $ahit ( @sorted_hits){
	print "map $trpt_id, $ahit->{sid}\n " if $_debug;
	if( $overlap = check_overlapping($trpt_id, $ahit->{sid})){
	    $best_hit = $ahit;
	    last;
	}
    }

    if($overlap){
	$qid2tidmapping{$trpt_id} = $best_hit;
    }else{
	warn ("Warn: $trpt_id, No hits found in vincinity\n");
	$qidnomapping{$trpt_id} = $sorted_hits[0];
    }
}



#====================================
# check if two transcript span overlap 
# in 1kb vincinity

sub check_overlapping{

    my $qID =shift;
    my $tID = shift;

    print "qID is $qID\ntID is $tID\n" if $_debug;

    my ($qchr, $qstart, $qend, $qstrand);
    if( $qID =~ /\|(\d+)\:(\d+)\:(\d+)\:([\d-]+)\:/){
        ($qchr, $qstart, $qend, $qstrand)   =($1, $2, $3, $4);
    }elsif($qID =~ /\|(\d+)\:(\d+)\:(\d+)\:([\d-]+)\|/){
	($qchr, $qstart, $qend, $qstrand)   =($1, $2, $3, $4);
    }else{
	warn("ERROR: bad query id: $qID\n");
	next;
    };

    my ($tchr, $tstart, $tend, $tstrand);
    if( $tID =~ /\|(\d+)\:(\d+)\:(\d+)\:([\d-]+)\|/){
        ($tchr, $tstart, $tend, $tstrand)   =($1, $2, $3, $4);
    }elsif($tID =~ /\|(\d+)\:(\d+)\:(\d+)\:([\d-]+)\:/ ){
	($tchr, $tstart, $tend, $tstrand)   =($1, $2, $3, $4);
    }else{
	warn("ERROR: bad target id: $tID\n");
	next;
    };

    if( $qchr ne $tchr or
	$qstrand ne $tstrand){
	return 0;
    }

    my $qinterval = new Number::Interval( Min => $qstart, IncMin => 1, 
					  Max => $qend, IncMax => 1 );
    my $tinterval = new Number::Interval( Min => $tstart, IncMin => 1, 
					  Max => $tend, IncMax => 1 );

    my $status = $qinterval->intersection ($tinterval);
    
    return $status;


}



#====================================
# print out the mapping
# blastpID format:
# query_trpt_stable_id|chr:chr_start:chr_end:chr_strand:translation_stable_id:gene_stable_id
# transcript_stable_id(v3) -> transcript_stable_id (v1)
# translation_stable_id(v3) -> translation_stable_id (v1)
# gene_stable_id(v3) -> gene_stable_id (v1)
#
# %qid2tidmapping;

sub print_id_mapping{

    my $index = 1;

    my %transcript_stable_id_mapping;
    my %translation_stable_id_mapping;
    my %gene_stable_id_mapping;

    foreach my $qID (keys %qid2tidmapping) {
	
	my $tID = $qid2tidmapping{$qID}->{sid};
	my $attrib = join "\t", 
	(map{ $qid2tidmapping{$qID}->{$keys->{$_}}}(3..14)); 

	my ($qtranscript_stable_id, $qtranslation_stable_id, $qgene_stable_id);
	if( $qID =~ /^([\w.]+)\|.+\:([\w.]+)\:([\w.]+)$/){
	    ($qtranscript_stable_id, $qtranslation_stable_id, $qgene_stable_id)
		=($1, $2, $3);
	}elsif($qID =~ /^([\w.]+)\|.+\|([\w.]+)\|([\w.]+)$/ ){
	    ($qtranscript_stable_id, $qtranslation_stable_id, $qgene_stable_id)
		=($1, $2, $3);
	}
	else{
	    warn("ERROR: bad query id: $qID\n");
	    next;
	};

	my ($ttranscript_stable_id, $ttranslation_stable_id, $tgene_stable_id);
	if( $tID =~ /^([\w.]+)\|.+\:([\w.]+)\:([\w.]+)$/){
	    ($ttranscript_stable_id, $ttranslation_stable_id, $tgene_stable_id)
		=($1, $2, $3);
	}elsif($tID =~ /^([\w.]+)\|.+\|([\w.]+)\|([\w.]+)$/ ){
	    ($ttranscript_stable_id, $ttranslation_stable_id, $tgene_stable_id)
		=($1, $2, $3);
	}else{
	    warn("ERROR: bad hit id: $tID\n");
	    next;
	};

	$transcript_stable_id_mapping{ $qtranscript_stable_id } = "$ttranscript_stable_id\t$attrib";
	$translation_stable_id_mapping{ $qtranslation_stable_id } = $ttranslation_stable_id;
	$gene_stable_id_mapping{ $qgene_stable_id } = $tgene_stable_id;
    }	

    my %outfile_sufix2table = (
	't' => \%transcript_stable_id_mapping,
	'tl' => \%translation_stable_id_mapping,
	'g' => \%gene_stable_id_mapping,
	);
    
    for my $suffix (keys %outfile_sufix2table){
	open my $ofh, '>', $outfile.".$suffix" or die "cannot open $outfile.$suffix to write";
	my $hashref = $outfile_sufix2table{$suffix};
	
	print $ofh join ("\n", (map{ $_."\t".$hashref->{$_}}  keys %{$hashref}) );
	print $ofh "\n";
	close $ofh;
	system("gzip -f $outfile.$suffix") and warn("Fail to gzip  $outfile.$suffix");
    }
    
}


sub print_id_nomapping{


    my $suffix = 'failed';
    open my $ofh, '>', $outfile.".$suffix" or die "cannot open $outfile.$suffix to write";
    
    my $i=2;
    map{ my $qid = $_;
	 my $besthit=$qidnomapping{$qid};
	 print $ofh join ("\t", ($qid, 
				 (map{ $besthit->{$keys->{$_}} } (2..14))
				 ));
	 print $ofh "\n";
	 
    } keys %qidnomapping;

    close $ofh;
    system("gzip -f $outfile.$suffix") and warn("Fail to gzip  $outfile.$suffix");
    
    
}



=head1 NAME

transcript_mapping_from_blastp_result.pl

=head1 SYNOPSIS

  transcript_mapping_from_blastp_result.pl    -infile input_file_name -outfile outfilename


Options:

  -h|--help        Show brief help and exit
  -v|--version     Show version and exit
  -infile          input file name, the one with blastp tab format 6, assuming the file is in the current working directory
  -outfile        output file name, write to current working directory

=head1 DESCRIPTION

To parse blastp result in tab format (6), and find the best target gene for each
query gene, first by sequence similarity, second by location vincinity.

=head1 SEE ALSO

perl.

=head1 AUTHOR

Sharon Wei E<lt>weix@cshl.eduE<gt>.

=head1 COPYRIGHT

Copyright (c) 2005 Cold Spring Harbor Laboratory

This library is free software;  you can redistribute it and/or modify 
it under the same terms as Perl itself.

=cut

