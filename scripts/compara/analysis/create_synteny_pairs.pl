#! /usr/bin/env perl

=head1 SYNOPSIS

create_synteny_pairs.pl  [options] 
 
	generate pairs of genomes we want build synteny from orthologs, dicot agaist arabidopsis and grape, monocot against rice and maize

 Options:
    --help              help message
    --man               full documentation
    --list         the file with speices list, each species a line, with dicot/monocot annotation


=cut


use Getopt::Long;
use Pod::Usage;
use Readonly;

Readonly::Hash my %syntney_refs => (
	dicot => [ qw(arabidopsis_thaliana vitis_vinifera) ], 
	monocot => [ qw(oryza_sativa zea_mays) ]
	);

my ($spfile,  $man, $help);

GetOptions (
    "list=s"          => \$spfile,
    "man"             => \$man,
    "help"            => \$help,
    ) or pod2usage(2);

pod2usage(2) if $man or $help;

my $fh;

open $fh, $spfile or die "Cannot open file $spfile to read";

my @species = map{ chomp; $_;} <$fh>;
my %species_2category;

map{ my ($k, $v)=split /\s+/; $species_2category{$k}=$v;} @species;

close $fh;

my @pairwise_pairs = pairup(\%species_2category);

#map { my $p=join "\t", @$_; print " pair is $p\n"} @pairwise_pairs;
#
for my $pair ( @pairwise_pairs ){

	my @args = @$pair;
	my $cmd = "$args[0]\t$args[1]";
	
        print "$cmd\n";
	#system($cmd) unless $test;

}

sub pairup{
	my $sp2cat = shift;
	my @pairs;

	for my $sp (keys %{$sp2cat}){
		my $cat = $sp2cat->{$sp};
		my $refs = $syntney_refs{$cat};

		foreach my $aref( @$refs ){
			push @pairs, [$sp, $aref];
		}
	}
	return @pairs;
}
