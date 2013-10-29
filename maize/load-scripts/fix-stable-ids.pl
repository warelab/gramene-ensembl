#!/usr/local/bin/perl

=head1 NAME

fix-stable-ids.pl - Resets all accession version numbers

=head1 SYNOPSIS

perl fix-stable-ids.pl [options]

Options:
 -h --help
 -m --man
 -r --registry_file
 -s --species
 -c --conf

=head1 OPTIONS

Resets the archiving scheme in the database. All accessions become versioned.

B<-h --help>
  Print a brief help message and exits.

B<-m --man>
  Print man page and exit

B<-r --registry>
  Use this Ensembl registry file for database connection info.
  Default is <ENSEMBLHOME>/conf/ensembl.registry

B<-s --species>
  Use this species entry from the registry file for the BAC database [REQUIRED].
  
B<-v --verbose>
  Provides verbose output.
  
=head1 DESCRIPTION

B<This program> 

  Resets the archiving scheme in the database. All accessions become versioned.
  
  Maintained by Shiran Pasternak <shiran@cshl.edu>

=cut

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

use Bio::EnsEMBL::Registry;

use English;
use FindBin qw($Bin);
use File::Basename qw( dirname );

use Log::Log4perl;

use vars qw($BASEDIR);

use Readonly;
use English qw( -no-match-vars );

use Memoize;
memoize('version_for');

my $logger = Log::Log4perl->get_logger('ensembl');

Readonly my $SUSPICIOUS_CLONES_SQL => <<SQL;
select substr(name, 1, 8) as name from seq_region where name like '%.'
SQL

Readonly my $REMOVE_SEQ_REGION_SQL => <<SQL;
delete from seq_region where name = ?
SQL

Readonly my $RENAME_SEQ_REGION_SQL => <<SQL;
update seq_region set name = ? where name = ?
SQL

Readonly my $UPDATE_SEQ_REGION_NAME_SQL => <<SQL;
update seq_region set name = ? where seq_region_id = ?
SQL

Readonly my $REMAP_CONTIGS_SQL => <<SQL;
update assembly
   set asm_seq_region_id = (select seq_region_id from seq_region
                             where name = ?)
 where cmp_seq_region_id in (select seq_region_id from seq_region
                              where name like ?)
SQL

BEGIN {

    # Set the perl libraries
    $BASEDIR = dirname($Bin);

    # unshift @INC, $BASEDIR . '/ensembl-live/ensembl/modules';
    Log::Log4perl::init($BASEDIR . '/conf/log4perl.conf');
}

MAIN: {
    my $help = 0;
    my $man  = 0;
    my ($species_name, $registry_file, $verbose);

    my $db_adaptor = undef;

    GetOptions(
        "help|?"     => \$help,
        "man"        => \$man,
        "species=s"  => \$species_name,
        "registry=s" => \$registry_file,
        "verbose"    => \$verbose,
    ) or pod2usage(2);

    pod2usage(-verbose => 2) if $man;
    pod2usage(1) if $help;

    # Validate file paths
    $registry_file ||= $BASEDIR . '/conf/SiteDefs.pm';

    map {
        -e $_ || (warn("File $_ does not exist\n")    && pod2usage(1));
        -r $_ || (warn("Cannot read $_\n")            && pod2usage(1));
        -f $_ || (warn("File $_ is not plain-text\n") && pod2usage(1));
        -s $_ || (warn("File $_ is empty\n")          && pod2usage(1));
    } $registry_file;

    $species_name || (warn("Need a --species\n") && pod2usage(1));
    Bio::EnsEMBL::Registry->load_all($registry_file);

    $db_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($species_name, 'core')
        || (warn("No core DB for $species_name set in $registry_file: $@\n")
        && pod2usage(1));

    fix_gene_names($db_adaptor);
}

=pod

=head2 fix_gene_names
    Ensures that all genes have the same versioned_accession as the prefix

=cut
sub fix_gene_names {
    my ($db_adaptor) = @_;
    
    my $gene_adaptor = $db_adaptor->get_GeneAdaptor();
    
}


=pod

=head2 update_slice_name
    Updates the slice name

=cut

sub update_slice_name {
    my ($adaptor, $slice, $name) = @_;

    update($adaptor->dbc->db_handle, $UPDATE_SEQ_REGION_NAME_SQL, $name,
        $slice->get_seq_region_id);
}

=pod

=head2 version_for
    Returns the version of this clone

=cut

sub version_for {
    my ($clone) = @_;

    return $clone->get_all_Attributes('acc-version')->[0]->value();
}

=pod

=head2 get_all_clones
    Returns all clones

=cut

sub get_all_clones {
    my ($adaptor) = @_;

    my $slice_adaptor = $adaptor->get_SliceAdaptor();

    my $clones = $slice_adaptor->fetch_all('clone');
    $logger->debug("Got @{[ scalar @$clones ]} clones");
    return $clones;
}

=pod

=head2 suspicious_clones
    Returns clone names that have bad archives

=cut

sub suspicious_clones {
    my ($dbh) = @_;

    return @{ $dbh->selectcol_arrayref($SUSPICIOUS_CLONES_SQL) };
}

=pod

=head2 remove_empty_clone
    Removes clone that doesn't have any contigs associated with it

=cut

sub remove_empty_clone {
    my ($dbh, $clone) = @_;

    # print "\tRemoving with SQL => $REMOVE_SEQ_REGION_SQL ($clone)\n";

    update($dbh, $REMOVE_SEQ_REGION_SQL, $clone);
}

=pod

=head2 rename_faulty_clone_name
    Renames the clone to be appropriate

=cut

sub rename_faulty_clone_name {
    my ($dbh, $clone) = @_;

    # print "\tRename clone $clone. => $clone\n";

    update($dbh, $RENAME_SEQ_REGION_SQL, $clone, "$clone.");
}

=pod

=head2 remap_contigs_to_fixed_clone
    Remap contigs from this clone in the assembly table

=cut

sub remap_contigs_to_fixed_clone {
    my ($dbh, $clone) = @_;

    update($dbh, $REMAP_CONTIGS_SQL, $clone, "$clone-Contig\%");
}

=pod

=head2 update
    Updates the DB

=cut

sub update {
    my ($dbh, $sql, @params) = @_;

    my $statement = $dbh->prepare($sql);
    $statement->execute(@params);
}
