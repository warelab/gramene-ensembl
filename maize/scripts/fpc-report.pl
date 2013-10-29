#!/usr/local/bin/perl

=head1 NAME

fpc-report.pl - Reports on sequenced clones on the FPC map

=head1 SYNOPSIS

perl fpc-report.pl [options]

Options:
 -a --(no)all
 -h --help
 -m --man
 -r --registry_file
 -c --conf

=head1 OPTIONS

B<-a --(no)all>
  List all clones, sequenced or not. B<--noall> suppresses other clones (DEFAULT: true).

B<-h --help>
  Print a brief help message and exits.

B<-m --man>
  Print man page and exit

B<-r --registry>
  Use this Ensembl registry file for database connection info.
  Default is <ENSEMBLHOME>/conf/ensembl.registry

B<-v --verbose>
  Provides verbose output.
  
=head1 DESCRIPTION

B<This program> 

  Prints out a tab-delimited report of sequenced clones on the FPC map
  
  Maintained by Shiran Pasternak <shiran@cshl.edu>

=cut

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

use Bio::EnsEMBL::Registry;

use FindBin qw($Bin);
use File::Basename qw( dirname );

use Log::Log4perl;

use vars qw($BASEDIR);

use Readonly;
use English qw( -no-match-vars );

use Perl6::Slurp;
use File::Spec;

my $logger = Log::Log4perl->get_logger('ensembl');

Readonly my $SQLDIR => "$BASEDIR/sql";

BEGIN {

    # Set the perl libraries
    $BASEDIR = dirname($Bin);    
    
    # unshift @INC, $BASEDIR . '/ensembl-live/ensembl/modules';
    Log::Log4perl::init($BASEDIR . '/conf/log4perl.conf');
}

MAIN: {
    my $help = 0;
    my $man  = 0;
    my $all  = 0;
    my ($species_name, $registry_file, $verbose, $all_clones);

    my $db_adaptor = undef;

    GetOptions(
        "help|?"     => \$help,
        "man"        => \$man,
        "registry=s" => \$registry_file,
        "all|a!"     => \$all_clones
    ) or pod2usage(2);
    
    $all_clones = 1 unless defined $all_clones;
    $species_name = 'Zea_mays';

    pod2usage(-verbose => 2) if $man;
    pod2usage(1) if $help;

    # Validate file paths
    $registry_file ||= "$BASEDIR/conf/ensembl-staging.registry";

    map {
        -e $_ || (warn("File $_ does not exist\n")    && pod2usage(1));
        -r $_ || (warn("Cannot read $_\n")            && pod2usage(1));
        -f $_ || (warn("File $_ is not plain-text\n") && pod2usage(1));
        -s $_ || (warn("File $_ is empty\n")          && pod2usage(1));
    } $registry_file;

    $species_name || (warn("Need a species\n") && pod2usage(1));
    Bio::EnsEMBL::Registry->load_all($registry_file);

    $db_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($species_name, 'core')
        || (warn("No core DB for $species_name set in $registry_file: $@\n")
        && pod2usage(1));

    print_report($db_adaptor, $all_clones);
}

=pod

=head2 print_report
    Prints the report

=cut
sub print_report {
    my ($adaptor, $all_clones) = @_;
    
    my $sqlfile = 'fpc_sequenced_bacs.sql';
    $sqlfile = 'all_fpc_clones.sql' if $all_clones;
    
    my $print_row = sub { print join("\t", @_), "\n"; };
    
    my $fullpath = File::Spec->catdir($SQLDIR, $sqlfile);
    my $sql = slurp($fullpath);

    my $statement = $adaptor->dbc->db_handle->prepare($sql);
    $statement->execute();
    $print_row->(@{ $statement->{NAME_lc} });
    while (my @row = $statement->fetchrow_array) {
        $print_row->(@row);
    }

}
