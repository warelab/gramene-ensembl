#!/usr/local/bin/perl -w 

=head1 NAME

load_assembly_from_fpc.pl - Populates a core Ensembl DB with FPC data

=head1 SYNOPSIS

perl load_assembly_from_fpc.pl [options] fpc_file

Options:

 -h --help
 -m --man
 -r --registry_file
 -s --species
 -c --cb_to_bp
 -n --no_insert

=head1 OPTIONS

Reads the B<fpc_file>, and uses its data to populate A. the
seq_region table and B. the assembly table of the core Ensembl DB.

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

B<-n --no_insert>
  Do not make changes to the database. Useful for debug.


=head1 DESCRIPTION

B<This program> 

  Populates a core Ensembl DB with FPC assembly data.
  Could do with additional description!
  
  A good place to look for clone sequences is;
  /usr/local/data/fpc/<species>/<species>.fpc
  
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
  #unshift @INC, $BASEDIR.'/bioperl-live'; # Need a recent version
}

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::CoordSystem;
use Bio::EnsEMBL::Slice;
use Bio::MapIO;

use vars qw( $ENS_DBA $FPC_MAP $SCALING_FACTOR $I );
BEGIN{  #Argument Processing
  my( $help, $man, $no_insert );
  
  my( $species, $file, $cb_to_bp );
  GetOptions
      (
       "help|?"          => \$help,
       "man"             => \$man,
       "no_insert"       => \$no_insert,
       "species=s"       => \$species,
       "registry_file=s" => \$file,
       "cb_to_bp=s"      => \$cb_to_bp,
       ) or pod2usage(2);
  pod2usage(-verbose => 2) if $man;
  pod2usage(1) if $help;

  # CB to BP scaling foctor. Calculated from AGI site. 
  # Can we get this programatically?
  $SCALING_FACTOR = $cb_to_bp || 4096;

  $file    ||= $BASEDIR.'/conf/ensembl.registry';
  my $fpc_file = shift @ARGV;
  $fpc_file || ( warn( "Need the path to an FPC file\n" ) && pod2usage(1));

  map{
    -e $_ || ( warn( "File $_ does not exist\n" )    && pod2usage(1) );
    -r $_ || ( warn( "Cannot read $_\n" )            && pod2usage(1) );
    -f $_ || ( warn( "File $_ is not plain-text\n" ) && pod2usage(1) );
    -s $_ || ( warn( "File $_ is empty\n" )          && pod2usage(1) );
  } $file, $fpc_file;
  
  $I = $no_insert ? 0 : 1;
  unless( $I ){ warn( "[INFO] Evaluation run - no db changes will be made\n")}

  # Load the ensembl file
  $species || ( warn( "Need a --species\n" ) && pod2usage(1) );
  Bio::EnsEMBL::Registry->load_all( $file );
  $ENS_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'core' );
  $ENS_DBA || ( warn( "No core DB for $species set in $file\n" ) &&
                pod2usage(1) );

  # Load the FPC file
  warn( "Found an FPC file: $fpc_file. Loading...\n" );
  my $mapio = new Bio::MapIO(-format  => "fpc",
                             -file    => "$fpc_file",
                             -readcor => 0,
                             -verbose => 0);
  $FPC_MAP = $mapio->next_map(); # Single map per FPC file
}

my( $chr_cs, $fpc_cs ) = &get_coord_systems();
my $sl_adapt = $ENS_DBA->get_adaptor('Slice');

# Loop through each chr and contig
foreach my $chr( &get_chromosome_names ){

  my $chr_length = &get_chromosome_length($chr);
  warn( "Processing chromosome $chr ($chr_length bp)\n" );

  my $chr_slice = &new_slice( $chr_cs,$chr,1,$chr_length );
  my $chr_dbid = 'XX';
  if( $I ){ $chr_dbid  = $sl_adapt->store($chr_slice) };

  my( $chr_start, $chr_end ) =  (1,1);
  foreach my $ctg( &get_contig_names($chr) ){

    my $ctg_length = &get_contig_length($ctg);
    $chr_end = $chr_start + $ctg_length -1;
    my $ctg_start = 1;

    my $ctg_name = "ctg${ctg}";
    warn( "  Processing $ctg_name - $ctg_length bp, ".
          "Chr$chr:$chr_start-$chr_end\n" );

    my $ctg_slice = &new_slice( $fpc_cs,$ctg_name,1,$ctg_length );
    my $ctg_dbid = 'XX';
    if( $I ){ $ctg_dbid  = $sl_adapt->store($ctg_slice) }
    
    &new_assembly( $chr_dbid,
                   $ctg_dbid,
                   $chr_start,
                   $chr_end,
                   $ctg_start,
                   $ctg_length );

    $chr_start = $chr_end + 1;

  }
}

#======================================================================
# Internal method to process FPC data into lists of clones grouped by
# chromosomes
use vars qw( $CHR_CTGS $CHR_LENGTHS $CTG_LENGTHS );
sub _process_fpc_data{
  $CHR_CTGS    ||= {};
  $CHR_LENGTHS ||= {};
  $CTG_LENGTHS ||= {};
  foreach my $ctgname( $FPC_MAP->each_contigid ){
    $ctgname || next; # Skip unnamed contigs
    my $ctgobj   = $FPC_MAP->get_contigobj($ctgname);

    my $ctgpos   = $ctgobj->position     || 0;
    my $ctgstart = $ctgobj->range->start || 0;
    my $ctgend   = $ctgobj->range->end   || 0;
    my $ctgbasepair = ( $ctgend - $ctgstart ) * $SCALING_FACTOR;
    $ctgpos += 0; # Force into numerical context
    my $chr = $ctgobj->group || '';
    $chr ||= 'UNKNOWN';
    #unless( $chr =~ /^\d+/ ){ $chr="UNKNOWN" }

    $CHR_CTGS->{$chr} ||= [];
    $CHR_LENGTHS->{$chr} ||= 0;

    push @{$CHR_CTGS->{$chr}}, [$ctgname, $ctgpos];
    $CHR_LENGTHS->{$chr} += $ctgbasepair;
    $CTG_LENGTHS->{$ctgname} = $ctgbasepair;
  }
  foreach my $chr( keys %$CHR_CTGS ){
    $CHR_CTGS->{$chr} = [ map{ $_->[0] }
                          sort{ ($a->[1]<=>$b->[1]) || 
                                ($a->[0]<=>$b->[0]) } @{$CHR_CTGS->{$chr}} ]
  }
  return $CHR_CTGS;
}

#----------------------------------------------------------------------
# Returns a list of all chromosomes in the FPC file
sub get_chromosome_names{
  $CHR_CTGS || &_process_fpc_data;
  return sort keys( %{$CHR_CTGS} );
}

#----------------------------------------------------------------------
# Returns the length in bp of a given chromosome
sub get_chromosome_length{
  $CHR_LENGTHS || &_process_fpc_data;
  my $chr = shift || die( "Need a chromosome name" );
  return $CHR_LENGTHS->{$chr};
}

#----------------------------------------------------------------------
# Returns a list of all contigs in the FPC file for a given chromosome,
# sorted in order of chromosome location
sub get_contig_names{
  my $chr = shift || die( "Need a chromosome name" );
  $CHR_CTGS || &_process_fpc_data;
  return @{ $CHR_CTGS->{$chr} || []};
}

#----------------------------------------------------------------------
# Returns the length in bp of a given contig
sub get_contig_length{
  $CTG_LENGTHS || &_process_fpc_data;
  my $ctg = shift || die( "Need a contig name" );
  return $CTG_LENGTHS->{$ctg};
}

#----------------------------------------------------------------------
# Creates a new slice based on the given cs_name, name, start and end
sub new_slice{

  return Bio::EnsEMBL::Slice->new
      (
       -coord_system      => shift || die( "Need a coord system" ),
       -seq_region_name   => shift || die( "Need a seq_region_name" ),
       -start             => shift || die( "Need a start coord" ),
       -end               => shift || die( "Need an end coord" ),
       -strand            => 1,
       -version           => 1,
       #-seq_region_length => $length,
       #-seq               => $seq->seq(), # Passed directly to store method
       );
}

#----------------------------------------------------------------------
# Inserts a new record into the assembly table based on the args
sub new_assembly{
  my $asm_seq_region_id = shift || die( "need an asm_seq_region_id" );
  my $cmp_seq_region_id = shift || die( "need acmp_seq_region_id" );
  my $asm_start = shift || die( "need an asm_start" );
  my $asm_end   = shift || die( "need an asm_end" );
  my $cmp_start = shift || die( "need an cmp_start" );
  my $cmp_end   = shift || die( "need an cmp_end" );
  my $ori = shift || 1;

  my $sql = qq"
INSERT INTO assembly
       (asm_seq_region_id, cmp_seq_region_id,
        asm_start, asm_end, cmp_start, cmp_end, ori) 
VALUES (?,?,?,?,?,?,?)";

  my $sth = $ENS_DBA->dbc->prepare($sql);
  my $rv;
  if( $I ){
    $rv = $sth->execute( $asm_seq_region_id, 
                         $cmp_seq_region_id,
                         $asm_start, $asm_end, 
                         $cmp_start, $cmp_end, $ori ) 
        || die( $sth->errstr );
  }
  return $rv;

}

#----------------------------------------------------------------------
# Retrieves the chromosome and fpc coord system adaptors
# Creates them if they do not exist
sub get_coord_systems{
  my $csa = $ENS_DBA->get_adaptor('CoordSystem');
  my $rank = 0;
  my %coord_systems;
  my @cs_names = ( 'chromosome','fpc' );
  foreach my $cs_name( @cs_names ){
    $coord_systems{$cs_name} = $csa->fetch_by_name($cs_name);
    unless( $coord_systems{$cs_name} ){
      $rank = ++$rank;
      $coord_systems{$cs_name} = Bio::EnsEMBL::CoordSystem->new
          (
           -name           => $cs_name,
           -rank           => $rank,
           -default        => 1,
           #$cs_name eq 'chromosome' ? () : (-sequence_level=>1),
           );
      print( "Storing coord_system: $cs_name\n" );
      if( $I ){
        $csa->store($coord_systems{$cs_name});
      }
    }
    $rank = $coord_systems{$cs_name}->rank;
  }

  # Populate meta table assembly.mapping entry (chromosome|fpc)
  my $meta = $ENS_DBA->get_MetaContainer();
  my $found;
  foreach my $map( @{$meta->list_value_by_key('assembly.mapping')} ){
    my @bits = split( /\|/, $map );
    if( $bits[0] = $cs_names[0] and $bits[1] eq $cs_names[1] ){
      $found ++;
    }
  }
  unless( $found ){
    my $sql = "
INSERT INTO meta( meta_key, meta_value )
VALUES ( ?, ? )";
    my $sth = $ENS_DBA->dbc->prepare($sql);
    if( $I ){
      my $rv = $sth->execute( 'assembly.mapping', join( '|', @cs_names ) )
          || die( $sth->errstr );
    }
  }
  
  return( map{$coord_systems{$_}} @cs_names );
}
