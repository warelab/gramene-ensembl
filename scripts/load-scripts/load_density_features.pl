#!/usr/local/bin/perl -w

=head1 NAME

load_density_features.pl - Populates a core Ensembl DB with density features

=head1 SYNOPSIS

perl load_density_features.pl [options]

Options:
 [B<-h>|B<--help>]
 [B<-m>|B<--man>]
 [B<-r>|B<--registry_file> I<file>]
 [B<-s>|B<--species> I<species>]
 [B<-l>|B<--logic_name> I<name>]
 [B<-v>|B<--verbose>]
 [B<-d>|B<--debug>]
 [B<-n>|B<--no_insert>]

=head1 OPTIONS

B<-h|--help>
  Print a brief help message and exits.

B<-m|--man>
  Print man page and exit

B<-r|--registry_file> I<file>
  Use this Ensembl registry file for database connection info.
  Default is <ENSEMBLHOME>/conf/ensembl.registry

B<-s|--species> I<species>
  Use this species entry from the registry file [REQUIRED].

B<-l|--logic_name> I<name>
  Indicates the analysis logic_name of the features to calculate
  densities from. In the case of misc_features, the misc_set.code
  should be used instead. [REQUIRED]

B<-v|--verbose>
  Print verbose output

B<-d|--debug>
  Print very verbose output

B<-n|--no_insert>
  Do not insert anything into the database. I.e. read-only mode.
  Useful for debug.

=head1 DESCRIPTION

B<This program> 

  Populates a core Ensembl DB with density features. Density features
  represent the density of a source feature type (logic_name) within
  given seq_region (e.g. chromosome) bins. Density feature data are
  used for the 'density' tracks displayed on mapview ideogram images.

B<The Ensembl Registry>

  The database connection details for the Ensembl database must be
  configured in an ensembl registry file.

  An example registry file would be;
  ---
  use Bio::EnsEMBL::DBSQL::DBAdaptor; 
  use Bio::EnsEMBL::Registry;
  # The Ensembl core database
  Bio::EnsEMBL::DBSQL::DBAdaptor->new
  ( '-species' => 'rice', 
    '-group'   => 'core', 
    '-dbname'  => 'oryza_sativa_core_28_17', );
  ---

  The script would then be called using
  shell> perl load_protein_features.pl -s=rice -l=whatever


B<Restoring the database>

  If the script bombs out half-way through, your database will be in a
  partial loaded state (i.e. in a bit of a mess). Here is some SQL to
  help put it right;

  TODO: Add the SQL!

Maintained by Will Spooner <whs@ebi.ac.uk>

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
  unshift @INC, $BASEDIR.'/ensembl-live/ensembl-compara/modules';
  unshift @INC, $BASEDIR.'/bioperl-live';
}

# Ensembl modules
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DensityType;
use Bio::EnsEMBL::DensityFeature;

use vars qw( $ENS_DBA $LOGIC_NAME $V $VV $INSERT );
BEGIN{  #Argument Processing
  my $help=0;
  my $man=0;
  my( $species, $file, $group, $no_insert );
  GetOptions
      ( 
        "help|?"            => \$help,
        "man"               => \$man,
        "species=s"         => \$species,
        "group=s"           => \$group,
        "registry_file=s"   => \$file,
        "logic_name=s"      => \$LOGIC_NAME,
        "no_insert"         => \$no_insert,
        "verbose"           => \$V,
        "debug"             => \$VV,
        )
      or pod2usage();
  pod2usage(-verbose => 2) if $man;
  pod2usage(1) if $help;

  $INSERT= $no_insert ? 0 : 1; # Put stuff in the database?

  # Load the ensembl file
  $LOGIC_NAME || ( print( "Need a --logic_name\n" ) && pod2usage() );
  $species || ( print( "Need a --species\n" ) && pod2usage() );
  $file    ||= $BASEDIR.'/conf/ensembl.registry';
  -e $file || ( print( "File $file does not exist\n" )    && pod2usage() );
  -r $file || ( print( "Cannot read $file\n" )            && pod2usage() );
  -f $file || ( print( "File $file is not plain-text\n" ) && pod2usage() );
  -s $file || ( print( "File $file is empty\n" )          && pod2usage() );

  Bio::EnsEMBL::Registry->load_all( $file );
  $ENS_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'core' );
  $ENS_DBA || ( print( "No core DB for $species set in $file\n" ) &&
                pod2usage() );
}


#
# Look up the logic_name, and determine corresponding feature type
#
my @feature_types;
my %block_counts;

# Check analysis.logic_name
my $aa = $ENS_DBA->get_adaptor("Analysis"); 
foreach my $type( $aa->feature_classes() ){

  #print "get feature class $type\n";
  #Skip Affy_Feature because the reported error "Table 'oryza_sativa_core_39_22.affy_feature' doesn't exist at /usr/local/gramene_ensembl-HEAD/ensembl-live/ensembl/modules/Bio/EnsEMBL/DBSQL/AnalysisAdaptor.pm"

  next if $type eq 'AffyFeature';
  foreach my $analysis( @{$aa->fetch_all_by_feature_class($type)} ){
    if( uc($analysis->logic_name) eq uc($LOGIC_NAME) ){ print "matched $LOGIC_NAME, $type\n";
      push( @feature_types, $type );

      my $table = join( '_', map lc, ( $type =~ /([A-Z][a-z]+)/g ) );
      my $count_sql = qq(
SELECT COUNT(*) 
FROM  $table 
WHERE analysis_id=? );
      
      my $sth = $ENS_DBA->dbc()->prepare( $count_sql );
      $sth->execute($analysis->dbID) || die( $sth->errstr );
      my $feature_count = $sth->fetchrow_array();
      print( "[INFO] $feature_count ${type}s found for $LOGIC_NAME\n" ) if $V;
      $block_counts{$type} = $feature_count >> 1; # Halve?
      
      last;
    }
  }
}

# Check misc_set.code
my $msa = $ENS_DBA->get_adaptor("MiscSet");
if( my $set = $msa->fetch_by_code( $LOGIC_NAME ) ){
  my $type = "MiscFeature";
  push( @feature_types, $type );

  my $count_sql = qq(
SELECT COUNT(*)
FROM  misc_feature_misc_set
WHERE misc_set_id=? );

  my $sth = $ENS_DBA->dbc()->prepare( $count_sql );
  $sth->execute($set->dbID) || die( $sth->errstr );
  my $feature_count = $sth->fetchrow_array();
  print( "[INFO] $feature_count ${type}s found for $LOGIC_NAME\n" ) if $V;
  $block_counts{$type} = $feature_count >> 1; # Halve?                      
}

# Quit if logic name not found
@feature_types or
    ( print( "\n[*DIE] $LOGIC_NAME does not correspond to a valid ".
             "analysis.logic_name or misc_set.code\n\n" )
      && pod2usage());

# Shout if logic name maps to more than one feature type
@feature_types > 1 and
    ( print( "$LOGIC_NAME maps to multiple feature types: ".
             join(', ', @feature_types )."\n" ) );

# Take first hit, and display type if verbose
my $feature_type = $feature_types[0];
print( "Logic_name $LOGIC_NAME ".
       "maps to features of type $feature_type\n" ) if $V;

#
# Check that database holds an assembly
#
ASSM_CHECK:{
  my $sth = $ENS_DBA->dbc()->prepare( "select count(*)  from seq_region" );
  $sth->execute();
  my ( $seq_region_count ) = $sth->fetchrow_array();
  $seq_region_count > 0 or die( "No seq_regions for core database!\n" );
}

#
# Get the adaptors needed;
#
my $dfa           = $ENS_DBA->get_adaptor('DensityFeature');
my $dta           = $ENS_DBA->get_adaptor('DensityType');
my $slice_adaptor = $ENS_DBA->get_adaptor('Slice');

#
# block size estimation
#

my $top_slices = $slice_adaptor->fetch_all('chromosome');
my $genome_size;
for my $slice ( @$top_slices ) {
  $genome_size += $slice->length();
}
my $block_count = $block_counts{$feature_type};
my $block_size = int( $genome_size / $block_count );
	

#
# Create new analysis entry based lon logic_name
#
my $analysis = new Bio::EnsEMBL::Analysis 
    (-program     => "logic_name_calc.pl",
     -database    => "ensembl",
     -gff_source  => "logic_name_calc.pl",
     -gff_feature => "density",
     -logic_name  => "Density_$LOGIC_NAME");
$aa->store( $analysis ) if $INSERT;

#
# Now the actual feature calculation
#
# Create new density type.
#
$analysis = $aa->fetch_by_logic_name("Density_$LOGIC_NAME") 
    if $INSERT;
my $dt = Bio::EnsEMBL::DensityType->new(-analysis   => $analysis,
                                        -block_size => $block_size,
                                        -value_type => 'sum');
if ($INSERT){

# Delete the old computational records for this analysis
# in density_feature and density_type tables 
    my $density_types = $dta->fetch_all_by_logic_name("Density_$LOGIC_NAME");
    foreach my $density_type(@$density_types){
	my $dt_id = $density_type->dbID();
	$dfa->db->dbc->do("delete from  density_feature where density_type_id = $dt_id");
	$dta->db->dbc->do("delete from  density_type where density_type_id = $dt_id");
    }
    $dta->store($dt) ;
}

my ( $current_start, $current_end );

my $class = join('', map{ucfirst($_)} split('_', $feature_type) );
my $getter_method = "get_all_${class}s";

foreach my $slice (@$top_slices){

  $current_start = 1;

  my @density_features=();

  print( "$LOGIC_NAME densities for ".$slice->seq_region_name().
         " with block size $block_size\n" ) if $V;

  while($current_start <= $slice->end()) {
    $current_end = $current_start+$block_size-1;
    if( $current_end > $slice->end() ) {
      $current_end = $slice->end();
    }
    
    my $sub_slice = $slice->sub_Slice( $current_start, $current_end );
    
    my $count = scalar( @{$sub_slice->$getter_method($LOGIC_NAME)} );
    push @density_features, Bio::EnsEMBL::DensityFeature->new
        (-seq_region    => $slice,
	 -start         => $current_start,
	 -end           => $current_end,
	 -density_type  => $dt,
	 -density_value => $count);
    
    $current_start = $current_end + 1;
  }
  $dfa->store(@density_features) if $INSERT;
  print( "Created ". scalar @density_features. 
         " $LOGIC_NAME density features.\n" ) if $V;
  print_features(\@density_features) if $VV;
}


#
# helper to draw an ascii representation of the density features
#
sub print_features {
  my $features = shift;

  return if(!@$features);

  my $sum = 0;
  my $length = 0;
#  my $type = $features->[0]->{'density_type'}->value_type();

  print("\n");
  my $max=0;
  foreach my $f (@$features) {
    if($f->density_value() > $max){
      $max=$f->density_value();
    }
  }
  if( !$max ) { $max = 1 };

  foreach my $f (@$features) {
    my $i=1;
    for(; $i< ($f->density_value()/$max)*40; $i++){
      print "*";
    }
    for(my $j=$i;$j<40;$j++){
      print " ";
    }
    print "  ".$f->density_value()."\t".$f->start()."\n";
  }
}




  


