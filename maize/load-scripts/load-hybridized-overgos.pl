#!/usr/local/bin/perl

=head1 NAME

load-hybridized-overgos.pl - Loads overgo-clone associations into a core Ensembl DB
Adding misc_set manually for misc_feature type

=head1 SYNOPSIS

perl load-hybridized-overgos.pl [options] fpc.file

Options:
 -h --help
 -m --man
 -r --registry_file
 -s --species

=head1 OPTIONS

Reads the B<fpc.file>, and uses its clones to load the misc_attrib table.

B<-h --help>
  Print a brief help message and exits.

B<-m --man>
  Print man page and exit

B<-r --registry_file>
  Use this Ensembl registry file for database connection info.
  Default is <ENSEMBLHOME>/conf/ensembl.registry

B<-s --species>
  Use this species entry from the registry file [REQUIRED].

=head1 DESCRIPTION

B<This program> 

  Loads overgo-clone associations into a core Ensembl DB
  
  Maintained by Shiran Pasternak <shiranp@cshl.edu>

=cut

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper qw(Dumper);    # For debug

use English;
use Carp;
use DBI;
use FindBin qw($Bin);
use File::Basename qw( dirname );

use vars qw($BASEDIR);

BEGIN {

    # Set the perl libraries
    $BASEDIR = dirname($Bin);
    unshift @INC, $BASEDIR . '/ensembl-live/ensembl/modules';
}

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::MiscFeature;
use Bio::EnsEMBL::Attribute;
use Bio::EnsEMBL::SimpleFeature;
use Bio::EnsEMBL::MiscSet;
use Bio::MapIO;

use Bio::EnsEMBL::Map::Marker;
use Bio::EnsEMBL::Map::MarkerSynonym;
use Bio::EnsEMBL::Map::MarkerFeature;

use vars qw($ENS_DBA $FPC_MAP);

my $help = 0;
my $man  = 0;
my ($species_name, $registry_file);
GetOptions(
	   "help|?"          => \$help,
	   "man"             => \$man,
	   "species=s"       => \$species_name,
	   "registry_file=s" => \$registry_file,
	   )
    or pod2usage(2);
pod2usage(-verbose => 2) if $man;
pod2usage(1) if $help;

# Validate file paths
$registry_file ||= $BASEDIR . '/conf/SiteDefs.pm';
my $fpc_file = shift @ARGV;
$fpc_file
    || (warn("Need the path to an overgo hits file\n") && pod2usage(1));

map {
    -e $_ || (warn("File $_ does not exist\n")    && pod2usage(1));
    -r $_ || (warn("Cannot read $_\n")            && pod2usage(1));
    -f $_ || (warn("File $_ is not plain-text\n") && pod2usage(1));
    -s $_ || (warn("File $_ is empty\n")          && pod2usage(1));
} $registry_file, $fpc_file;

# Load the ensembl file
$species_name || (warn("Need a --species\n") && pod2usage(1));
Bio::EnsEMBL::Registry->load_all($registry_file);
$ENS_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $species_name, 'core' );
$ENS_DBA    || (warn("No core DB for $species_name set in $registry_file\n")
        && pod2usage(1));

# Load the FPC file
warn( "Found an FPC file: $fpc_file. Loading...\n" );
my $mapio = new Bio::MapIO(-format  => "fpc",
			   -file    => "$fpc_file",
			   -readcor => 0,
			   -verbose => 0);
$FPC_MAP = $mapio->next_map(); # Single map per FPC file

# species info
my $meta    = $ENS_DBA->get_MetaContainer();
my $species = $meta->get_Species
    || die("Cannot find the species in the meta table of the DB");
my $common_name = $species->common_name
    || die(
	   "Cannot find the species common name in the meta table of the DB");

$common_name = ucfirst($common_name);

###########
# Prepare some Ensembl adaptors
my $misc_feature_adaptor = $ENS_DBA->get_adaptor('MiscFeature');
my $attribute_adaptor    = $ENS_DBA->get_adaptor('Attribute');

# store associations for each marker in the FPC datastruc
my @markers = $FPC_MAP -> each_markerid();
warn( "  Num markers: ",scalar(@markers) );

my $associations = 0;
foreach my $marker (@markers){
#    last if ($associations > 50);
    my $markerobj = $FPC_MAP -> get_markerobj($marker);
    my @clones = $markerobj -> each_cloneid();
    my @contigs = $markerobj -> each_contigid();
    
#    print "$marker\t";
#    foreach my $contig (@contigs){
#	print "$contig;";
#    }
#    print "\t";
#    foreach my $clone (@clones){
#	print "$clone;";
#    }
#    print "\n";
    
    my $marker_attribute = get_overgo_attribute($marker);
    for my $clone_name (@clones){
	my $clone_features
	    = $misc_feature_adaptor->fetch_all_by_attribute_type_value('name', $clone_name);

	if (scalar @$clone_features != 1) {
	    warn
		"Unexpected feats for <$clone_name>; expected 1, got @{[scalar @$clone_features]}";
	    next;
	}
	my $clone_feature = $clone_features->[0];
	my $ctg = $clone_feature->get_all_Attributes('superctg');	
	my $ctg_name = ($ctg->[0]) -> value();
	
	if (clone_already_hybridized_with_overgo($clone_feature, $marker_attribute)){
	    warn "Clone $clone_name already associated with $marker\n";
	} 
	else {
	    warn
		"Associating $marker with: @{[$clone_feature->get_scalar_attribute('name')]}\n";
#	    $clone_feature->add_Attribute($marker_attribute);
#	    $attribute_adaptor->store_on_MiscFeature($clone_feature,
#						     [$marker_attribute]);
	    $associations++;
	}
    }
}


#======================================================================

=pod

=head2 get_overgo_attribute
    Ensure the existence of an overgo attribute type

=cut

sub get_overgo_attribute {
    my ($attribute_value) = @_;

    # my ($dba) = @_;

    # my $attribute_adaptor = $dba->get_adaptor('Attribute');
    return Bio::EnsEMBL::Attribute->new(
        -CODE        => 'clone_marker',
        -NAME        => 'Clone Marker',
        -DESCRIPTION => 'A marker with which a clone is associated',
        -VALUE       => $attribute_value,
    );
}

=pod

=head2 clone_already_hybridized_with_overgo
    Determines if the attribute already exists

=cut

sub clone_already_hybridized_with_overgo {
    my ($clone_feature, $overgo_attribute) = @_;
    my @attributes =
        @{$clone_feature->get_all_Attributes($overgo_attribute->code)};
    return (grep { $overgo_attribute->value() eq $_->value() } @attributes) > 0;
}

#======================================================================
1;


