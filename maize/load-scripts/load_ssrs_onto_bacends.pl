#!/usr/local/bin/perl -w

=head1 NAME

load_ssrs_onto_bacends.pl - Annotates bacs in a core Ensembl DB with SSRs

=head1 SYNOPSIS

perl load_ssrs_onto_bacends.pl [options] data_file

Options:
 [B<-h>|B<--help>]
 [B<-m>|B<--man>]
 [B<-r>|B<--registry_file> I<file>]
 [B<-v>|B<--verbose>]
 [B<-d>|B<--debug>]
 [B<-n>|B<--no_insert>]

=head1 OPTIONS

Reads the B<data_file>, and uses its mapping data to assign SSRs
(consensus + length) to clones in the misc_features table (identified
by misc_set.name == "BAC map"). Association is via the misc_attrib
table, with an attrib_type.code of "bacend_ssr".

B<-h|--help>
  Print a brief help message and exits.

B<-m|--man>
  Print man page and exit

B<-r|--registry_file> I<file>
  Use this Ensembl registry file for database connection info.
  Default is <ENSEMBLHOME>/conf/ensembl.registry

B<-s|--species>
  Use this species entry from the registry file [REQUIRED].

B<-v|--verbose>
  Print verbose output

B<-d|--debug>
  Print very verbose output

B<-n|--no_insert>
  Do not insert anything into the database. I.e. read-only mode.
  Useful for debug.


=head1 DESCRIPTION

B<This program> 

B<The data_file>

  This file contains the mappings between the GenBank accession for
  the BACend with the SSR data. The columns are;
   - bacend_accession (e.g. CL000917),
   - ssr_no (e.g. 1),
   - ssr_type(di,tri,tetra) (e.g. 3),
   - ssr_consensus (e.g. ttg)
   - no_of_repeats (e.g. 12),
   - ssr_start (e.g. 513),
   - ssr_end (e.g. 536),
   - bacend_length (e.g. 1065)


  A good place to look for this file is
  /usr/local/data/analysis/<species>/<species>_ssr_to_bacend.dat
  
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
}

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::MiscFeature;
use Bio::EnsEMBL::MiscSet;
use Bio::EnsEMBL::Attribute;
use Bio::EnsEMBL::Analysis;

use vars qw( $ENS_DBA $DATA_FILE $V $VV $INSERT );
BEGIN{  #Argument Processing
  my $help=0;
  my $man=0;
  my( $species, $regfile, $group, $no_insert );
  GetOptions
      ( 
        "help|?"            => \$help,
        "man"               => \$man,
        "registry_file=s"   => \$regfile,
        "species=s"         => \$species,
        "no_insert"         => \$no_insert,
        "verbose"           => \$V,
        "debug"             => \$VV,
        )
      or pod2usage();
  pod2usage(-verbose => 2) if $man;
  pod2usage(1) if $help;
  $VV && $V++; # Use verbose if debug

  $DATA_FILE = shift || 
      ( print( "\n[*DIE] Need a data_file\n\n" ) && pod2usage() );

  $species || 
      ( print( "\n[*DIE] Need a species\n\n" ) && pod2usage() );

  $regfile    ||= $BASEDIR.'/conf/ensembl.registry'; # Def registry file
  
  $INSERT= $no_insert ? 0 : 1; # Put stuff in the database?

  # Validate files
  foreach my $file( $DATA_FILE, $regfile ){
    -e $file || ( print( "\n[*DIE] File $file does not exist\n\n" ) 
                  && pod2usage() );
    -r $file || ( print( "\n[*DIE] Cannot read $file\n\n" ) 
                  && pod2usage() );
    -f $file || ( print( "\n[*DIE] File $file is not plain-text\n\n" )
                  && pod2usage() );
    -s $file || ( print( "\n[*DIE] File $file is empty\n\n" )
                  && pod2usage() );
  }
  print( "[INFO] Loading data from $DATA_FILE\n" ) if $V;

  Bio::EnsEMBL::Registry->load_all( $regfile );
  $ENS_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'core' );
  $ENS_DBA || ( print( "\n[*DIE] No core DB for $species in $regfile\n\n" ) &&
                pod2usage() );
}

# Open the DATA_FILE
open( FH, $DATA_FILE );

# Prepare adaptors
my $mfa = $ENS_DBA->get_adaptor('MiscFeature');
my $aa  = $ENS_DBA->get_adaptor('Attribute');
# Read each record
my $found = 0;
my $unfound = 0;
my %bacend_cache = ();
my $header = <FH>;
while( my $line = <FH> ){
  my( $bacend_name, 
      $ssr_no,
      $ssr_type,
      $ssr_consensus,
      $num_repeats,
      $ssr_start,
      $ssr_end,
      $bacend_length ) = split( /\s+/, $line );

  my $feature = $bacend_cache{$bacend_name};

  unless( defined( $feature ) ){
    print "[INFO] Processing BACend $bacend_name\n" if $VV;

    # Find the misc feature corresponding to the bac name. Clone names
    # can be prefixed with library. Try both with and without
    my @feats = @{$mfa->fetch_all_by_attribute_type_value('bacend',
                                                          $bacend_name)};
    
    # Take the first feature that has a misc_set of bac_map
    foreach my $feat( @feats ){
      foreach my $set( @{$feat->get_all_MiscSets('bac_map')} ){
        $feature = $feat;
        last;
      }
      $feature && last;
    }
    
    if( $feature ){
      # Found BAC
      print( "[INFO] BACend $bacend_name found: @feats\n" ) if $VV;
      $bacend_cache{$bacend_name} = $feature;
      $found ++;
    } else {
      # Could not find BAC!
      print( "[WARN] BACend $bacend_name not found in DB\n" ) if $VV;
      $bacend_cache{$bacend_name} = 0;
      $unfound ++;
    }    
  }
  $feature || next;

  my $repeat = "($ssr_consensus)$num_repeats";

  my $attrib = Bio::EnsEMBL::Attribute->new(-VALUE => $repeat,
                                            -CODE  => 'bacend_ssr',
                                            -NAME  => 'BACend SSR');

  $aa->store_on_MiscFeature( $feature, [$attrib] ) if $INSERT;

}

print( "[INFO] processed ".($found+$unfound).
       " BACs. $found found, $unfound missing\n" );

#======================================================================
1;
