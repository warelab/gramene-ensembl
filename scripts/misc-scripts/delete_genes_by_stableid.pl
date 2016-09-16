#!/usr/bin/env perl
=head1 NAME

delete_genes_by_stableid.pl - delete genes and their transcript and associated attribs.

	

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

delete_genes_by_stableid.pl  [options] 
 
 Options:
    --help		help message
    --man		full documentation
    --registry          the registry file for database connections
    --speceis 		which species to dump
    --list	        list of gene stable Ids whose gene are to be removed


=head1 OPTIONS

=over 4


    
=item B<--registry>

The registry file for ensembl databases

=item B<--species> 

supply the species name whose genes/transcripts are to be deleted

=item B<--list> 

File of  list of gene stable Ids whose gene are to be removed

=item B<--help> 

print a help message and exit

=item B<--man> 

print documentation and exit

    
=back



=cut

my ($species, $registry, $lfile, $write);


{  #Argument Processing
  my $help=0;
  my $man=0;

  GetOptions( "help|?"=>\$help,"man"=>\$man
	      ,"write" => \$write
	      ,"list=s"=>\$lfile
	      ,"species=s"=>\$species
	      ,"registry=s"=>\$registry

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
my $slice_adaptor = $ENS_DBA->get_SliceAdaptor; 
#my $assembly_type=$ENS_DBA->get_MetaContainer()->get_default_assembly();

my $genes;

if ($lfile){
    open my $fh, $lfile or die "Cannot open gene stable id file for reading";
    push @$genes, map{ print; chomp;$gene_adaptor->fetch_by_stable_id($_) } <$fh>;
    close $fh;
}

my %count;


foreach my $gene (@$genes) {
    my $sid = $gene->stable_id;
  print "! gene_stable_id=$sid\n";
  
  $count{genes2delete}++;
 
  if($write){
      eval {$gene_adaptor->remove($gene)};
      if( $@ ){
	  print STDERR "ERROR: cannot removed gene $sid, $@";
      }
  }
  
  
}

for my $k (sort keys %count){
    print "$k = $count{$k}\n";
}




