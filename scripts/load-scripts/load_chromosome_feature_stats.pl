#!/usr/local/bin/perl -w

=head1 NAME

load_chromosome_feature_stats.pl - Populates the seq_region_attribute
table of an Ensembl DB with counts of all features that it knows
about for each top-level sequence region.

=head1 SYNOPSIS

perl load_chromosome_feature_stats.pl [options]

Options:
 -h --help
 -m --man
 -r --registry_file
 -s --species
 -l --logic_name
 -v --verbose
 -d --debug
 -n --no_insert

=head1 OPTIONS

B<-h --help>
  Print a brief help message and exits.

B<-m --man>
  Print man page and exit

B<-r --registry_file>
  Use this Ensembl registry file for database connection info.
  Default is <ENSEMBLHOME>/conf/ensembl.registry

B<-s --species>
  Use this species entry from the registry file [REQUIRED].

B<-l --logic_name> 
  Restrict to a given analysis logic_name(s) or misc_set class. 
  If omitted, all applicable analyses and misc_sets are processed.

B<-v --verbose>
  Print verbose output

B<-d --debug>
  Print very verbose output

B<-n --no_insert>
  Do not insert anything into the database. I.e. read-only mode.
  Useful for debug.

=head1 DESCRIPTION

B<This program> 

  Could do with more description

Maintained by Will Spooner <whs@ebi.ac.uk>

=head1 WARNINGS

Should probably delete all counts from seq region feature when 
recalculating all counts.

=head1 SEE ALSO

The Ensembl scripts for generating per-chromosome lengths,
gene counts and so on:
ensembl/misc-scripts/density_feature/seq_region_stats.pl

There is also a Gramene SOP here;
http://gwiki.gramene.org/Ensembl_Core_Seq_Region_Stats

=cut

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper qw(Dumper); # For debug
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
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Attribute;

use vars qw( @LOGIC_NAMES $ENS_DBA $V $VV $INSERT );
BEGIN{  #Argument Processing
  my $help=0;
  my $man=0;
  my( $species, $file, $group, $no_insert );
  GetOptions
      ( 
        "help|?"            => \$help,
        "man"               => \$man,
        "species=s"         => \$species,
        "registry_file=s"   => \$file,
        "logic_name=s"      => \@LOGIC_NAMES,
        "no_insert"         => \$no_insert,
        "verbose"           => \$V,
        "debug"             => \$VV,
        )
      or pod2usage(2);
  pod2usage(-verbose => 2) if $man;
  pod2usage(1) if $help;

  # Load the ensembl file
  $species || ( print( "Need a --species\n" ) && pod2usage(1) );
  $file    ||= $BASEDIR.'/conf/ensembl.registry';
  -e $file || ( print( "File $file does not exist\n" )    && pod2usage(1) );
  -r $file || ( print( "Cannot read $file\n" )            && pod2usage(1) );
  -f $file || ( print( "File $file is not plain-text\n" ) && pod2usage(1) );
  -s $file || ( print( "File $file is empty\n" )          && pod2usage(1) );
  Bio::EnsEMBL::Registry->load_all( $file );
  $ENS_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'core' );
  $ENS_DBA || ( print( "No core DB for $species set in $file\n" ) &&
                pod2usage() );
  
  $INSERT= $no_insert ? 0 : 1; # Put stuff in the database?

  print( "Loading seq_region_attributes for $species;\n" );
  foreach my $dba( $ENS_DBA ){
    print( "  Database: ".$dba->dbc->dbname ."@". $dba->dbc->host .
           ":". $dba->dbc->port ."\n");
  }
}

my $slice_adaptor    = $ENS_DBA->get_adaptor('Slice');
my $attrib_adaptor   = $ENS_DBA->get_adaptor('Attribute');
my $analysis_adaptor = $ENS_DBA->get_adaptor('Analysis');
my $miscset_adaptor  = $ENS_DBA->get_adaptor('MiscSet');

my $top_slices = $slice_adaptor->fetch_all( "chromosome" ) ||
  $slice_adaptor->fetch_all( "toplevel" );


my %analyses    = ();
my %logic_names = ();
foreach my $class( "DnaAlignFeature", 
                   "ProteinAlignFeature", 
                   "MarkerFeature", #MarkerFeature seem empty for now
                   "SimpleFeature",
                   "Gene",
                   "PredictionTranscript") {
  $analyses{$class} = $analysis_adaptor->fetch_all_by_feature_class($class);
  map{ $logic_names{ $_->logic_name} = $class } @{$analyses{$class}};
}
foreach my $class( "MiscFeature" ){
  $analyses{$class} = $miscset_adaptor->fetch_all;
  map{ $logic_names{ $_->code} = $class } @{$analyses{$class}};
}

# Restrict logic_names to those requested (if any)
if( @LOGIC_NAMES ){
  my %requested_logic_names = ();
  foreach my $lname( @LOGIC_NAMES ){
    my $class = $logic_names{$lname} ||
      (print( "\nLogic name $lname not found in database\n" ) && pod2usage() );
    $requested_logic_names{$lname} = $class;
  }
  %logic_names = %requested_logic_names;
}

#AttributeAdaptor doesn't have the appropriate remove method (only all
#	attrib of slice) So we do sql:
$attrib_adaptor->db->dbc->do(
    "delete from seq_region_attrib where attrib_type_id in
        (select attrib_type_id from attrib_type where code in ("
    . join(",", map { "'${_}Count'" } keys %logic_names)
    .") )" 
);

# Dump the data...
foreach my $slice (@$top_slices) {
  my $seq_region_name = $slice->seq_region_name;
  print "Processing seq_region $seq_region_name\n";
  my @attribs;
  foreach my $class( keys %analyses ){
    print("  Feature class $class\n");
    my $getter_method = "get_all_${class}s";

    foreach my $analysis( @{$analyses{$class}} ){
      my $lname;
      my $label;
      if( $analysis->isa('Bio::EnsEMBL::Analysis') ){
        $lname = $analysis->logic_name;
        $label = $analysis->display_label || $lname;
      } elsif ( $analysis->isa('Bio::EnsEMBL::MiscSet') ){
        $lname = $analysis->code;
        $label = $analysis->name || $lname;
      } else {
        die( "Cannot process objects of type ",ref($analysis) );
      }
      $logic_names{$lname} || next;
      my $desc  = "Total Number of $label features";
      my $count = scalar( @{$slice->$getter_method($lname)} );
      print("    $desc on $seq_region_name = $count\n" );
      push @attribs, Bio::EnsEMBL::Attribute->new
          (-NAME        => "$label Count",
           -CODE        => substr($lname.'Count', 0, 40), #a hack, the db limit is 40
           -VALUE       => $count,
           -DESCRIPTION => $desc);
    }
  }

  $attrib_adaptor->store_on_Slice($slice, \@attribs) if $INSERT;
#  print_chromo_stats([$slice]);
}



sub print_chromo_stats {
  my $chromosomes = shift;

  foreach my $chr (@$chromosomes) {
    print "\nchromosome: ",$chr->seq_region_name(),"\n";
    foreach my $attrib (@{$chr->get_all_Attributes()}) {
      print "  ", $attrib->name(), ": ", $attrib->value(), "\n";
    }
  }
}


1;


