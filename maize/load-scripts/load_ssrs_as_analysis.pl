#!/usr/local/bin/perl -w

=head1 NAME

load_ssrs_as_analysis.pl - Populates a core Ensembl DB with ssr data 
    (as marker with analysis type)

=head1 SYNOPSIS

perl load_clones_from_fpc.pl [options] ssr_file

Options:
 -h --help
 -m --man
 -r --registry_file
 -s --species
 -c --cb_to_bp

=head1 OPTIONS

Reads the B<ssr_file>, and uses its bac end information to assign 
positions to ssrs in the marker_feature table. Should be run after
bac ends have been loaded to the misc_attrib table using 
load_bacends_onto_bacs.pl.

B<-h --help>
  Print a brief help message and exits.

B<-m --man>
  Print man page and exit

B<-r --registry_file>
  Use this Ensembl registry file for database connection info.
  Default is <ENSEMBLHOME>/conf/ensembl.registry

B<-s --species>
  Use this species entry from the registry file [REQUIRED].

B<-c --cb_to_bp>                                                                 
  Conversion factor to convert native (cb) fpc coordinated into basepairs.       
  Calculated by hand by comparison with coordinates shown on AGI browser.        
  Default: 4096.    

=head1 DESCRIPTION

B<This program> 

  Populates a core Ensembl DB with electronic ssrs
  
=cut

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper qw(Dumper); # For debug 

use DBI;
use FindBin qw( $Bin );
use File::Basename qw( dirname );

use vars qw( $BASEDIR );
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
use Bio::EnsEMBL::Analysis;

use Bio::EnsEMBL::Map::Marker;
use Bio::EnsEMBL::Map::MarkerSynonym;
use Bio::EnsEMBL::Map::MarkerFeature;

use vars qw( $ENS_DBA $SCALING_FACTOR );

# Arg Processing
my $help=0;
my $man=0;
my( $species_name, $file, $cb_to_bp );
GetOptions
    (
     "help|?"          => \$help,
     "man"             => \$man,
     "species=s"       => \$species_name,
     "registry_file=s" => \$file,
     "cb_to_bp=s"      => \$cb_to_bp,
     ) or pod2usage(2);
pod2usage(-verbose => 2) if $man;
pod2usage(1) if $help;

# Assign the scaling factor
$SCALING_FACTOR = $cb_to_bp;

# Validate file paths
$file    ||= $BASEDIR.'/conf/ensembl.registry';
my $ssr_file = shift @ARGV;
$ssr_file   || ( warn( "Need the path to an ssr file\n" ) && pod2usage(1));

map{
    -e $_ || ( warn( "File $_ does not exist\n" )    && pod2usage(1) );
    -r $_ || ( warn( "Cannot read $_\n" )            && pod2usage(1) );
    -f $_ || ( warn( "File $_ is not plain-text\n" ) && pod2usage(1) );
    -s $_ || ( warn( "File $_ is empty\n" )          && pod2usage(1) );
} $file, $ssr_file;

# Load the ensembl file
$species_name || ( warn( "Need a --species\n" ) && pod2usage(1) );
Bio::EnsEMBL::Registry->load_all( $file );
$ENS_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $species_name, 'core' );
$ENS_DBA || ( warn( "No core DB for $species_name set in $file\n" ) &&
	      pod2usage(1) );

my $meta    = $ENS_DBA->get_MetaContainer();
my $species = $meta->get_Species || 
    die( "Cannot find the species in the meta table of the DB" );
my $common_name = $species->common_name ||
    die( "Cannot find the species common name in the meta table of the DB" );

$common_name = ucfirst( $common_name );

###########
# Prepare some Ensembl adaptors
my $marker_analysis        = &fetch_analysis( "ssr_marker" );
my $marker_feature_adaptor = $ENS_DBA->get_adaptor('MarkerFeature');
#my $marker_adaptor         = $ENS_DBA->get_adaptor('Marker');
my $misc_feature_adaptor   = $ENS_DBA->get_adaptor('MiscFeature');
my $slice_adaptor          = $ENS_DBA->get_adaptor('Slice');

# Loop through each ssr
open (SSR, "$ssr_file");
while (my $line = <SSR>){
    chomp $line;
    
    my( $ssr_name,
	$bacend_name,
	$motif,
	$motif_length,
	$motif_count,
	$ssr_length,
	$ssr_start,
	$ssr_end,
	$bacend_length ) = split( /\t/, $line );

#    next unless ($ssr_name eq "cshr00056");
    
    # store name in the marker_syn table and add entry in the marker table
    my $ensmarker = Bio::EnsEMBL::Map::Marker->new();
    $ensmarker->display_MarkerSynonym
	( Bio::EnsEMBL::Map::MarkerSynonym->new(undef, "ssr_marker", $ssr_name) );
    $ensmarker->priority(100); #Ensure marker displays
    
    # get all features corresponding to the bac end (very slow)
    # check to see if there is ever more than one feature -- prob not.
    my @features = @{$misc_feature_adaptor->
			 fetch_all_by_attribute_type_value('bacend', $bacend_name)};
    
    # check to see if this bacend exists in the db, if not, report and bail
    my $feature = @features;
    if ($feature < 1){
	print STDERR "$ssr_name $bacend_name not found\n";
	next;
    }
    
    # get the start and end positions of the bac and the ctg
    my $inner_start;
    my $inner_end;
    my $ctg;
    foreach my $feature (@features){
	my @inner_start = @{$feature -> get_all_attribute_values('inner_start')};
	$inner_start    = $inner_start[0];
	my @inner_end   = @{$feature -> get_all_attribute_values('inner_end')};
	$inner_end      = $inner_end[0];
	my @ctg         = @{$feature -> get_all_attribute_values('superctg')}; 
	$ctg            = $ctg[0];
    }
    
    # estimate a marker location as average of start and end of the bac
    my $maploc = (($inner_start + $inner_end) / 2);
    
    # Create a slice from the FPC contig
    my $ctg_slice = $slice_adaptor->fetch_by_region(undef,$ctg);
    $ctg_slice || die( "FPC $ctg was not found in the Ensembl DB" );
    
    # Create the marker feature
    my $marker_feature = Bio::EnsEMBL::Map::MarkerFeature->new
        (
         undef, #DBid
         undef, #MarkerFeatureAdaptor
         $maploc,   # Start
         $maploc + $SCALING_FACTOR, # End - assume length=SCALING_FACTOR
         $ctg_slice,         # Slice
         $marker_analysis,   # Analysis
         undef,              # Marker ID?
         1,                  # map weight
         $ensmarker,         # Ensembl Marker
         );

    # transform to chrom coords for drawing code
    $marker_feature = $marker_feature->transform('chromosome') ||
        ( warn("No chr mapping for $ssr_name!") && next );
    
    print STDERR "Storing $ssr_name $bacend_name\n";
    
    $marker_feature_adaptor->store( $marker_feature );

}


# subs

# Returns an Analysis object; either fetched from the DB if one exists,
# or created fresh, in which case it is stored to the DB.
sub fetch_analysis{
    my $logic_name = shift || die("Need a logic_name" );
    my $db_file    = shift || '';
    my $adaptor = $ENS_DBA->get_adaptor('Analysis'); 
    
    my %args = ( -logic_name=>$logic_name, 
		 $db_file ? (-db_file=>$db_file) : () );
    
    my $analysis;
    if( $analysis = $adaptor->fetch_by_logic_name($args{-logic_name}) ){
	# Analysis found in database already; use this.
	return $analysis;
    }
    
    # No analysis - create one from scratch
    $analysis = Bio::EnsEMBL::Analysis->new(%args);
    $adaptor->store($analysis);
    return $analysis;
}

1;
