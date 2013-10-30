#!/usr/local/bin/perl 

=head1 NAME

cleanup_overkap_fgenesh - when fgenesh genes overlap each other, keep only the longest one
	

=cut


#use

#fetch_all_by_logic_name('foobar');

BEGIN {
    $ENV{'GrameneDir'} = '/usr/local/gramene/'; 
    $ENV{'GrameneEnsemblDir'} = '/usr/local/ensembl-live/'; 
}

# The first shall be last...
use lib map { $ENV{'GrameneDir'}."/$_" } qw { lib/perl };

use lib map { $ENV{'GrameneEnsemblDir'}."/$_" } 
        qw { bioperl-live modules ensembl/modules conf
	     ensembl-external/modules ensembl-draw/modules
	     ensembl-compara/modules };

use strict;
use warnings;
use List::MoreUtils qw{ distinct };


use Bio::SeqIO;
use Bio::EnsEMBL::DBLoader;
use Bio::EnsEMBL::Registry;
use Getopt::Long;
use Pod::Usage;


=head1 SYNOPSIS

dump_gene_coord.pl  [options] 
 
 Options:
    --help		help message
    --man		full documentation
    --registry          the registry file for database connections
    --speceis 		which species to dump
    --logic_name        logic name of the fgenesh analysis

=head1 OPTIONS

=over 4


=item B<--registry>

The registry file for ensembl databases

=item B<--species> 

supply the species name whose transcripts are to be dumped

=item B<--logic_name> 

 logic name of the fgenesh analysis

=item B<--help> 

print a help message and exit

=item B<--man> 

print documentation and exit

    
=back

=head1 ARGUMENTS

Input files with Gene Stable ids 

=cut

my ($species, $registry, $logic_name);

{  #Argument Processing
  my $help=0;
  my $man=0;

  GetOptions( "help|?"=>\$help,
	      "man"=>\$man,
	      "species=s"=>\$species,
	      "registry=s"=>\$registry,
	      "logic_name=s"=>\$logic_name,
	    )
    or pod2usage(2);

  pod2usage(-verbose => 2) if $man;
  pod2usage(1) if $help;

}


# Load the ensembl file
Bio::EnsEMBL::Registry->load_all( $registry );
my $ENS_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'core' );
$ENS_DBA || pod2usage( "\nNo core DB for $species set in $registry\n" );

my $exon_adaptor=$ENS_DBA->get_ExonAdaptor;
my $gene_adaptor=$ENS_DBA->get_GeneAdaptor;

my $cnt=0;
foreach my $gene (@{$gene_adaptor->fetch_all_by_logic_name($logic_name)}) {
  
#print "! $geneid\n";

    my $geneid = $gene->dbID;
    my $genename = $gene->display_id;
    my $seq_name = $gene->seqname;
    my $strand   = $gene->strand;
    my $start    = $gene->start;
    my $end      = $gene->end;

    my @genelist = grep {$_->dbID != $geneid} @{$gene->get_overlapping_Genes()};
    
    if( scalar @genelist >0  ){
	$cnt++;
	print "$cnt: ";
	print join '|', ($geneid, $genename, $seq_name, $strand, $start, $end);
	print ", overlapping other genes: ";
	print join ',', map{ join "|", ($_->dbID, $_->display_id)} @genelist;
	print "\n";
	    
    }

 }



=head1 AUTHOR

   Sharon Wei <weix@cshl.edu>

   Gramene Project (www.gramene.org)
   Ware Lab
   Cold Spring Harbor Laboratory

   This is free software and may be distributed on the same terms as Perl itself.

=cut

