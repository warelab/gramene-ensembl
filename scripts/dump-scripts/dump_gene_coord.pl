#!/usr/local/bin/perl 

=head1 NAME

dump_gene_coord.pl - dump coordinates of the genes in the input file
	

=cut

BEGIN {
    $ENV{'GrameneDir'} = '/usr/local/gramene/'; 
    $ENV{'GrameneEnsemblDir'} = '/usr/local/ensembl-live/'; 
}

# The first shall be last...
use lib map { $ENV{'GrameneDir'}."/$_" } qw ( lib/perl );
#use lib map { $ENV{'GrameneEnsemblDir'}."/ensembl-live/$_" } 
#        qw ( bioperl-live modules ensembl/modules ensembl-external/modules
#             ensembl-draw/modules ensembl-compara/modules );
use lib map { $ENV{'GrameneEnsemblDir'}."/$_" } 
        qw ( bioperl-live modules ensembl/modules conf
	     ensembl-external/modules ensembl-draw/modules
	     ensembl-compara/modules );

use strict;
use warnings;


#use Gramene::Config;
use Bio::SeqIO;
use Bio::EnsEMBL::DBLoader;
use Bio::EnsEMBL::Registry;
#use Bio::EnsEMBL::Gene;
use Getopt::Long;
use Pod::Usage;


=head1 SYNOPSIS

dump_gene_coord.pl  [options] 
 
 Options:
    --help		help message
    --man		full documentation
    --registry          the registry file for database connections
    --speceis 		which species to dump
    --index             geneid index seperated by ','

=head1 OPTIONS

=over 4


=item B<--registry>

The registry file for ensembl databases

=item B<--species> 

supply the species name whose transcripts are to be dumped

=item B<--index> 

index in the input file for gene names 0 based index for example "0,1,2"

=item B<--help> 

print a help message and exit

=item B<--man> 

print documentation and exit

    
=back

=head1 ARGUMENTS

Input files with Gene Stable ids 

=cut

my ($species, $registry, $index);

{  #Argument Processing
  my $help=0;
  my $man=0;

  GetOptions( "help|?"=>\$help,"man"=>\$man

	      ,"species=s"=>\$species
	      ,"registry=s"=>\$registry
	      ,"index=s"=>\$index
	    )
    or pod2usage(2);

  pod2usage(-verbose => 2) if $man;
  pod2usage(1) if $help;

}


# Load the ensembl file
Bio::EnsEMBL::Registry->load_all( $registry );
my $ENS_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'core' );
$ENS_DBA || pod2usage( "\nNo core DB for $species set in $registry\n" );

my $gene_adaptor=$ENS_DBA->get_GeneAdaptor;

foreach my $file (@ARGV) {
  #print "! $geneid\n";

  open my $fh, $file or die "Cannot open $file to read";
  
  my @rows = <$fh>;
  chomp @rows;
  close $fh;

  my $outfile = "$file.coord";
  open my $ofh, '>', $outfile or die "couldn't open $outfile to write";

  for my $line (@rows){

    my @fields   = split /\s+/, $line;
    #print "fields=@fields\n";
    my @geneids  = @fields[split ',',$index];
    #print "geneids=@geneids\n";
    my ($min_start, $max_end, @coords) = (0,0, undef);
    for my $g (@geneids){
      print "$g\n";
      my $gene;
      eval { $gene=$gene_adaptor->fetch_by_stable_id($g) };
      print STDERR "$@\n" and next if $@ || !$gene;
      
      my ($seq_region_name, $gs, $ge, $gstrand) = ($gene->seq_region_name,
						  $gene->seq_region_start,
						  $gene->seq_region_end,
						  $gene->seq_region_strand,);

      push @coords, ($g, $seq_region_name, $gs, $ge, $gstrand);

      if( $min_start == 0){
	$min_start = $gs < $ge ? $gs : $ge;
	$max_end   = $gs < $ge ? $ge : $gs;

      }else{

      for my $p ( $gs, $ge){
	
	$min_start = $p if $min_start > $p;
	$max_end   = $p if $max_end < $p;
      }
    }


      
    }
    print $ofh join ("\t", ($line, @coords, $min_start, $max_end));
    print $ofh "\n";
  }
  close $ofh;

}










=head1 AUTHOR

   Sharon Wei <weix@cshl.edu>

   Gramene Project (www.gramene.org)
   Ware Lab
   Cold Spring Harbor Laboratory

   This is free software and may be distributed on the same terms as Perl itself.

=cut

