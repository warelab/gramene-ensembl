#! /usr/bin/env perl

=head1 SYNOPSIS

submit_pariwise_pairs.pl  [options] 
 
	submit jobs for each pair of genome combinations, such as dumping orthologs

 Options:
    --help              help message
    --man               full documentation
    --list         the file with speices list, each species a line 
    --script       the submission script, usually customized for each job and store in working directory
    --test        print out submission command without actually submit


=cut


use Getopt::Long;
use Pod::Usage;
use Algorithm::Combinatorics qw( combinations );

my ($spfile, $script, $test, $man, $help);

GetOptions (
    "list=s"          => \$spfile,
    "script=s"        => \$script,
    "test"            => \$test,
    "man"             => \$man,
    "help"            => \$help,
    ) or pod2usage(2);

pod2usage(2) if $man or $help;

my $fh;

open $fh, $spfile or die "Cannot open file $spfile to read";

my @species = map{ chomp; $_;} <$fh>;

#map{print "got $_\n";} @species;

close $fh;

my @pairwise_pairs = combinations(\@species, 2);

#map { my $p=join "\t", @$_; print " pair is $p\n"} @pairwise_pairs;
#
for my $pair ( @pairwise_pairs ){

	my @args = @$pair;
	my $cmd = "qsub $script $args[0] $args[1] ";
	
        print "$cmd\n";
	system($cmd) unless $test;

}
