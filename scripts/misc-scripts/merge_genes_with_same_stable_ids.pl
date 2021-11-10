#!/usr/local/bin/perl 

=head1 NAME

merge_genes_with_same_ids.pl

=cut

BEGIN {
    $ENV{'GrameneDir'} = '/usr/local/gramene/'; 
    $ENV{'GrameneEnsemblDir'} = '/usr/local/ensembl-live/'; 
}

# The first shall be last...
use lib map { $ENV{'GrameneDir'}."/$_" } qw ( lib/perl );

use lib map { $ENV{'GrameneEnsemblDir'}."/$_" } 
        qw ( bioperl-live modules ensembl/modules conf
	     ensembl-external/modules ensembl-draw/modules
	     ensembl-compara/modules );

use strict;
use warnings;
use Scalar::Util qw (reftype);

#use Gramene::Config;
use Bio::SeqIO;
use Bio::EnsEMBL::DBLoader;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::TranscriptMapper;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Exon;

use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

=head1 SYNOPSIS

merge_genes_with_same_ids.pl  [options] 
 
 Options:
    --help		help message
    --man		full documentation
    --registry          the registry file for database connections
    --species 		which species to process
    --nowrite           test only , no change to the database
    --debug             debug mode

=head1 OPTIONS

=over 4


=item B<--registry>

    The registry file for ensembl databases

=item B<--species> 

    supply the species name whose translations are to be recalculated

=item B<--nowrite>

    test run without writing to database

=item B<--help> 

    print a help message and exit

=item B<--man> 

    print documentation and exit

=item B<--debug> 

   print out more debug information

=back

=head1 ARGUMENTS

=cut

my ($species, $registry);
my ($debug, $nowrite);
{  #Argument Processing
  my $help=0;
  my $man=0;

  GetOptions( "help|?"=>\$help,"man"=>\$man
	      ,"species=s"=>\$species
	      ,"registry=s"=>\$registry
	      ,"debug"=>\$debug
	      ,"nowrite"=>\$nowrite
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

print STDERR "#INFO DBbase connected is ", $ENS_DBA->dbc->dbname, "\n" if $debug;

my $dbh = $ENS_DBA->dbc->db_handle;

my $select_genes_by_biotype_sth = $dbh->prepare(qq{
  select stable_id, gene_id from gene where biotype=?
});

my $update_transcripts_sth =  $dbh->prepare(qq{
  update transcript set gene_id=? where gene_id=?
});

my $delete_gene_sth = $dbh->prepare("delete from gene where gene_id=?");

$select_genes_by_biotype_sth->execute("protein_coding");
my %protein_coding;
while (my $row = $select_genes_by_biotype_sth->fetchrow_arrayref) {
  my ($stable_id,$gene_id) = @$row;
  $protein_coding{$stable_id} = $gene_id;
}
$select_genes_by_biotype_sth->execute("non_coding");
my $updates=0;
while (my $row = $select_genes_by_biotype_sth->fetchrow_arrayref) {
  my ($stable_id,$gene_id) = @$row;
  if (exists $protein_coding{$stable_id}) {
    $updates++;
    print STDERR "merging transcripts of non_coding gene $stable_id ($gene_id)\n" if $debug;
    $update_transcripts_sth->execute($protein_coding{$stable_id}, $gene_id) unless $nowrite;
    print STDERR "deleting the non_coding gene $stable_id ($gene_id)\n" if $debug;
    $delete_gene_sth->execute($gene_id) unless $nowrite;
  }
}
$select_genes_by_biotype_sth->finish;
$update_transcripts_sth->finish;
$delete_gene_sth->finish;

$dbh->disconnect;

print STDERR "merged $updates genes\n";
  
########################## subroutines ######################################

__END__


=head1 OUTPUT


=head1 AUTHOR

   Andrew Olson <olson@cshl.edu>

   Gramene Project (www.gramene.org)
   Ware Lab
   Cold Spring Harbor Laboratory

   This is free software and may be distributed on the same terms as Perl itself.

=cut

