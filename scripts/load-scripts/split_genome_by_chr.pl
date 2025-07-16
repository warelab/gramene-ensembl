#!/lab/bin/perl 

=head1 NAME

tidy_fasta.pl - uniform the fasta files so that they can be compared by diff
                  the output files are fasta1.tidy, fasta2.tidy, the tidy process will
                  uniform all the charactor to capital and ordered by header 
                  alphalbetically

=cut


BEGIN {
    $ENV{'GrameneDir'} ||= '/usr/local/gramene/'; 
    $ENV{'GrameneEnsemblDir'} ||= '/usr/local/gramene_ensembl/'; 
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

#use Gramene::Config;
use DBI qw/:sql_types/;

use Bio::SeqIO;
use Bio::PrimarySeq;

use Bio::EnsEMBL::DBLoader;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::Registry;
use List::MoreUtils qw/each_array/;

=head1 SYNOPSIS

tidy_fasta.pl  [options] fasta1 fasta2 ...
 
 Options:
    --help		help message
    --man		full documentation
    --start_coord       the start coordinate
    -end_coord
    --seq_name          name of the sequence


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
  my ($seq_name, $outdir);

  GetOptions( 
	     "help|?"          => \$help,
	     "man"             => \$man,
	     "v|verbose"       => \$verbose,
	     "outdir=s"	       => \$outdir,
	    )
    or pod2usage(2);
  pod2usage(-verbose => 2) if $man;
  pod2usage(1) if $help;




$outdir ||= '.';

for my $f (@ARGV){
  split_genome_by_chr($f);
}


sub split_genome_by_chr{
  
  my ($fasta) = @_;
  my $fasta_seqio= Bio::SeqIO->new('-format' => 'fasta',
				   '-file' => $fasta)
    or die "can't open $fasta to read:$!";

  while( my $seq=$fasta_seqio->next_seq()){
    my $seqname = $seq->display_id;
    my $seqout= Bio::SeqIO->new( -format => 'Fasta', -file => ">$outdir/$seqname.fa");
    $seqout->write_seq( $seq );
  }
}

