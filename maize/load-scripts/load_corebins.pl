#!/usr/local/bin/perl -w

=head1 NAME

load_corebins.pl - Populates a core Ensembl DB with virtual core bin data

=head1 SYNOPSIS

perl load_corebins.pl [options] corebins.file

Options:
 -h --help
 -m --man
 -r --registry_file
 -s --species
 -n --no_insert

=head1 OPTIONS

Reads the B<corebins.file>, and uses its virtual bin data to 
load misc feature tables.

B<-h --help>
  Print a brief help message and exits.

B<-m --man>
  Print man page and exit

B<-r --registry_file>
  Use this Ensembl registry file for database connection info.

B<-s --species>
  Use this species entry from the registry file [REQUIRED].

B<-n --no_insert>
  Do not insert anything into the database. I.e. read-only mode.
  Useful for debug.

=head1 DESCRIPTION

B<This program> 

  Populates a core Ensembl DB with virtual core bins

=cut

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper qw(Dumper); # For debug 

use DBI;
use FindBin qw($Bin) ;
use File::Basename qw(dirname);

use vars qw($BASEDIR);
BEGIN{
  # Set the perl libraries
  $BASEDIR = dirname($Bin);
  unshift @INC, $BASEDIR.'/ensembl-live/ensembl/modules';
  unshift @INC, $BASEDIR.'/bioperl-live';
}

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::MiscFeature;
use Bio::EnsEMBL::MiscSet;
use Bio::EnsEMBL::Attribute;
use Bio::EnsEMBL::SimpleFeature;
#use Bio::EnsEMBL::Analysis;
#use Bio::MapIO;
#use Bio::SeqIO;

use vars qw($I $ENS_DBA);
my $help=0;
my $man=0;
my($species, $file, $no_insert);
GetOptions
    (
     "help|?"          => \$help,
     "man"             => \$man,
     "species=s"       => \$species,
     "registry_file=s" => \$file,
     "no_insert"       => \$no_insert,
     ) or pod2usage(2);
pod2usage(-verbose => 2) if $man;
pod2usage(1) if $help;

$I = $no_insert ? 0 : 1; # Put stuff in the database?

# Validate file paths
$file    ||= $BASEDIR.'/conf/ensembl.registry';
my $corebin_file = shift @ARGV;
$corebin_file   || (warn("Need the path to an FPC file\n") && pod2usage(1));

map{
    -e $_ || ( warn( "File $_ does not exist\n" )    && pod2usage(1) );
    -r $_ || ( warn( "Cannot read $_\n" )            && pod2usage(1) );
    -f $_ || ( warn( "File $_ is not plain-text\n" ) && pod2usage(1) );
    -s $_ || ( warn( "File $_ is empty\n" )          && pod2usage(1) );
} $file, $corebin_file;

# Load the ensembl file
$species || ( warn( "Need a --species\n" ) && pod2usage(1) );
Bio::EnsEMBL::Registry->load_all( $file );
$ENS_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'core' );
$ENS_DBA || ( warn( "No core DB for $species set in $file\n" ) &&
	      pod2usage(1) );

###########
# Prepare some Ensembl adaptors
my $sl_adapt = $ENS_DBA->get_adaptor('Slice');
my $f_adapt  = $ENS_DBA->get_adaptor('SimpleFeature');
my $mf_adapt = $ENS_DBA->get_adaptor('MiscFeature');

###########
# Establish the virtual core bin feature. Similar to analysis
my $misc_set_corebin = Bio::EnsEMBL::MiscSet->new
    (-code        => 'core_bins', # Must conform to Ensembl for display
     -name        => 'Virtual Core Bins',
     -description => 'IBM2/Neighbors and FPC intersection',
     -longest_feature => 10000000 );

##########                                                                         
# Process the core bin marker file and load vbins as misc features
open (VBINS, "$corebin_file");
while (my $line = <VBINS>){
    chomp $line;
    (my $bin, my $range) = split (/\t/, $line);
    (my $start, my $end) = split (/\-/, $range);
    (my $chr, my $binner) = split (/\./, $bin);

    my $chr_slice = $sl_adapt -> fetch_by_region ('chromosome', $chr);
    $chr_slice || die( "Chr $chr was not found in the Ensembl DB" );
    
    my $corebin_feature = Bio::EnsEMBL::MiscFeature -> new
	(
	 -start => $start,
	 -end   => $end,
	 -strand => 1,
	 -slice => $chr_slice,
	 );
    my $corebin_attrib = Bio::EnsEMBL::Attribute -> new
	(
	 -VALUE => $bin,
	 -CODE => 'name',
	 -NAME => 'name',
	 );
    
    $corebin_feature -> add_MiscSet ($misc_set_corebin);
    $corebin_feature -> add_Attribute ($corebin_attrib);
    $mf_adapt -> store ($corebin_feature) if $I;
}
close (VBINS);
