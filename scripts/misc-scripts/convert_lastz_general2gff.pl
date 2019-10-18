#!/env/bin perl 

=head1 NAME

	convert lastz out put file in general format to gff format
=cut


use strict;
use warnings;


use Getopt::Long;
use Pod::Usage;
use File::Glob;

use List::MoreUtils qw/each_array/;

=head1 SYNOPSIS

program  lastz_file_in_general_format
 
 Options:
    --help		help message
    --man		full documentation


=head1 OPTIONS

=over 4

    
=item B<--help> 

print a help message and exit

=item B<--man> 

print documentation and exit

   
=back

=head1 AUTHOR

Sharon Wei

=cut


    

{  #Argument Processing
  my $help=0;
  my $man=0;
  my $verbose;
  
  GetOptions( 
	     "help|?"          => \$help,
	     "man"             => \$man,
	     "v|verbose"       => \$verbose,

	    )
    or pod2usage(2);
  pod2usage(-verbose => 2) if $man;
  pod2usage(1) if $help;

}

my $input_file = shift;


open my $fh, $input_file or die "cannot open file $input_file";

#score  name1   strand1 size1   zstart1 end1    name2   strand2 size2   zstart2 end2    identity        idPct   coverage        covPct
#31878   1       +       308452471       1471239 1472085 1       +       80884392        6541    7387    583/846 68.9%   846/80884392    0.0%
#5420    1       +       308452471       9039092 9039430 1       +       80884392        15736   16074   192/338 56.8%   338/80884392    0.0%
#5541    1       +       308452471       9463801 9463891 1       +       80884392        19082   19172   76/90   84.4%   90/80884392     0.0%
#
#


my $score_idx = 0;
my $target_idx = 1;
my $target_start_idx = 4;
my $target_end_idx = 5;
my $query_start_idx = 9;
my $query_end_idx = 10;
my $source = 'lastz';
my $type = 'match';

while ( <$fh> ){

	next if ( /score/i );

	chomp;

	my @mapping = split ' ';

	my $target_chr = $mapping[$target_idx];
	my $target_start = $mapping[$target_start_idx];
	my $target_end = $mapping[$target_end_idx];
	my $score = $mapping[$score_idx];
	my $strand = $mapping[2] eq $mapping[7] ? '+' : '-';
	
	my %attribs;
	$attribs{t_size} = $mapping[3];
	$attribs{q_size} = $mapping[8];
	$attribs{q_chr} = $mapping[6];
	$attribs{q_start} = $mapping[9];
	$attribs{q_end} = $mapping[10];
	$attribs{identity} = $mapping[11];
	$attribs{id_pct} = $mapping[12];
	$attribs{coverage} = $mapping[13];
	$attribs{cov_pct} = $mapping[14];

	
	my $attribs = join ";", (map { "$_=$attribs{$_}"  } qw(q_chr q_size q_start q_end t_size identity id_pct coverage cov_pct));

	print join "\t", ( $target_chr, $source, $type, $target_start, $target_end, $score, $strand, '.', "$attribs\n");
}

close $fh;

