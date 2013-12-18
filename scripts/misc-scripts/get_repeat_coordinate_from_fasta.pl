#!/lab/bin/perl 

=head1 NAME

tidy_fasta.pl - uniform the fasta files so that they can be compared by diff
                  the output files are fasta1.tidy, fasta2.tidy, the tidy process will
                  uniform all the charactor to capital and ordered by header 
                  alphalbetically

=cut


BEGIN {
    $ENV{'GrameneDir'} ||= '/usr/local/gramene/'; 
    $ENV{'GrameneEnsemblDir'} ||= '/usr/local/ensembl-live/'; 
}

# The first shall be last...
use lib map { $ENV{'GrameneDir'}."/$_" } qw ( lib/perl );
use lib map { $ENV{'GrameneEnsemblDir'}."/ensembl-live/$_" } qw ( bioperl-live modules ensembl/modules ensembl-external/modules ensembl-draw/modules ensembl-compara/modules);
#use lib map { $ENV{'GrameneEnsemblDir'}."/$_" } qw ( bioperl-live modules ensembl/modules conf  ensembl-external/modules ensembl-draw/modules ensembl-compara/modules);



use strict;
use warnings;


use Getopt::Long;
use Pod::Usage;
use File::Glob;

use Gramene::Config;
use DBI qw/:sql_types/;

use Bio::SeqIO;
use Bio::PrimarySeq;

use Bio::EnsEMBL::DBLoader;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::Registry;
use List::MoreUtils qw/each_array/;

=head1 SYNOPSIS

get_repeat_coordinate_from_fasta.pl  [options] fasta1 fasta2 ...
 
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


    

  #Argument Processing
  my $help=0;
  my $man=0;
  my $verbose;
  my ($start, $end, $seq_name);

  GetOptions( 
	     "help|?"          => \$help,
	     "man"             => \$man,
	     "v|verbose"       => \$verbose,

	    )
    or pod2usage(2);
  pod2usage(-verbose => 2) if $man;
  pod2usage(1) if $help;




my $counter = 1 ;

for my $f (@ARGV){
    search_fasta($f);
}


sub search_fasta{
  
  my ($fasta) = @_;
  my $fasta_seqio= Bio::SeqIO->new('-format' => 'fasta',
				   '-file' => $fasta)
    or die "can't open $fasta to read:$!";

  my %header_bioseq_hash;
  while( my $seq=$fasta_seqio->next_seq()){
    
    my $seq_name = $seq->display_id;
    #print "seq_name=$seq_name\n";
    my $seq =  $seq->seq;
    while(  $seq =~ /[a-z]+/gc  ){
	my $end = pos($seq); #1-based
	my $start = $end - length($& )+ 1;
	print join "\t", ($counter++, $seq_name, $start, "$end\n");


    }

 #   	last;
  }

}

