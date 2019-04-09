#!/lab/bin/perl 

=head1 NAME

tidy_fasta.pl - uniform the fasta files so that they can be compared by diff
                  the output files are fasta1.tidy, fasta2.tidy, the tidy process will
                  uniform all the charactor to capital and ordered by header 
                  alphalbetically

=cut



use strict;
use warnings;


use Getopt::Long;
use Pod::Usage;
use File::Glob;

use DBI qw/:sql_types/;

use Bio::SeqIO;
use Bio::PrimarySeq;

use List::MoreUtils qw/each_array/;

=head1 SYNOPSIS

tidy_fasta.pl  [options] fasta1 fasta2 ...
 
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


for my $f (@ARGV){
  tidy_fasta($f);
}


sub tidy_fasta{
  
  my $fasta = shift;
  my $fasta_seqio= Bio::SeqIO->new('-format' => 'fasta',
				   '-file' => $fasta)
    or die "can't open $fasta to read:$!";

  my %header_bioseq_hash;
  while( my $seq=$fasta_seqio->next_seq()){
    
    my $header = uc($seq->display_id);
    my $seq_str = uc($seq->seq);
    $header_bioseq_hash{$header} = new Bio::Seq(
						-display_id => $header,
						-seq => $seq_str,
					       );
  }
  
  my $out= Bio::SeqIO->new('-format' => 'fasta',
		       '-file' => ">$fasta.tidy")
      or die "can't open $fasta.tidy for output:$!";
  
  for my $id(sort keys %header_bioseq_hash){
    $out->write_seq($header_bioseq_hash{$id});
  }
}

