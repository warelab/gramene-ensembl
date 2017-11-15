#!/usr/local/bin/perl -w

=head1 NAME

load_compara_synteny.pl - Populates a EnsEMBL compara database with synteny data from a text file

=head1 SYNOPSIS

perl load_compara_synteny.pl [options] synteny_file

Options:
 
 -h|--help
 -m|--man
 -r|--registry_file <file>
 -q|--query_species <registry_file_species>
 -t|--target_species <registry_file_species>
 -c|--compara <registry_file_compara_db_species>
 -v|--verbose
 -d|--debug
 -n|--no_insert

=head1 OPTIONS

B<-h|--help>
  Print a brief help message and exits.

B<-m|--man>
  Print man page and exit

B<-r|--registry_file> I<file>
  Use this Ensembl registry file for database connection info.
  Default is <ENSEMBLHOME>/conf/ensembl.registry

B<-q>|B<--query_species> I<registry_file_species>
  the species name in the registry_file for one of the synteny species, should correspond to the 1st species in the synteny file

B<-t>|B<--target_species> I<registry_file_species>
  the species name in the registry_file for the other synteny species, should correspond to the 2nd species in the synteny file

B<-c>|B<--compara> I<registry_file_species>
  the species name in the registry_file for the compara db

B<-v|--verbose>
  Print verbose output

B<-d|--debug>
  Print very verbose output

B<-n|--no_insert>
  Do not insert anything into the database. I.e. read-only mode.
  Useful for debug.

=head1 DESCRIPTION

B<This program> 

  Populates the compara database with precomputed syntenic blocks
  defined in file. At the moment the sepcies are hard-coded to  
  Zea_mays and Oryza_sativa_japonica (TODO: make this a command
  line arg). 

B<The synteny file>

  The format of the file is (whitespace separated);
sp1_seq_region_name sp1_start sp1_end sp2_seq_region_name sp2_start sp2_end

B<The Ensembl Registry>

  The database connection details for the Ensembl compara database
  must be configured in an ensembl registry file.

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

Created by Will Spooner <whs@ebi.ac.uk>
Modified by Sharon Wei <weix@cshl.edu>

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
  $BASEDIR = '/sonas-hs/ware/hpc/data/weix/ensembl_main/'; #'/usr/local'; #dirname($Bin);
  #$BASEDIR = dirname($Bin);
  unshift @INC, $BASEDIR.'/ensembl/modules';
  unshift @INC, $BASEDIR.'/ensembl-compara/modules';
  unshift @INC, $BASEDIR.'/bioperl-live';
}

# Ensembl modules
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DensityType;
use Bio::EnsEMBL::DensityFeature;

use vars qw( $ENS_DBA $SYNTENY_FILE $LOGIC_NAME $V $VV $INSERT 
             $BAND_TO_BASEPAIR);
my ($query_species, $target_species, $compara );

BEGIN{  #Argument Processing
  my $help=0;
  my $man=0;
  my( $species, $regfile, $group, $no_insert );
  GetOptions
      ( 
        "help|?"            => \$help,
        "man"               => \$man,
        "registry_file=s"   => \$regfile,
        "query_species=s"   => \$query_species,
        "target_species=s"  => \$target_species,
        "compara=s"         => \$compara,
        "no_insert"         => \$no_insert,
        "verbose"           => \$V,
        "debug"             => \$VV,
        )
      or pod2usage();
  pod2usage(-verbose => 2) if $man;
  pod2usage(1) if $help;
  $VV && $V++; # Use verbose if debug

  $SYNTENY_FILE = shift || ( print( "\n[*DIE] Need a SYNTENY_FILE\n\n" )
                             && pod2usage() );

  $regfile    ||= $BASEDIR.'/conf/ensembl.registry'; # Def registry file

  $query_species || (print "Need query species, the 1st species in the synteny file" && pod2usage() );
  $target_species || (print "Need target species, the 2nd species in the synteny file" && pod2usage() );
  $compara || (print "Need compara db species as is set in the registry file" && pod2usage() );

  $INSERT= $no_insert ? 0 : 1; # Put stuff in the database?

  $BAND_TO_BASEPAIR ||= 4096; # Unit conversion - band to basepair

  # Validate files
  foreach my $file( $SYNTENY_FILE, $regfile ){
    -e $file || ( print( "\n[*DIE] File $file does not exist\n\n" ) 
                  && pod2usage() );
    -r $file || ( print( "\n[*DIE] Cannot read $file\n\n" ) 
                  && pod2usage() );
    -f $file || ( print( "\n[*DIE] File $file is not plain-text\n\n" )
                  && pod2usage() );
    -s $file || ( print( "\n[*DIE] File $file is empty\n\n" )
                  && pod2usage() );
  }
  print( "[INFO] Reading synteny from $SYNTENY_FILE\n" ) if $V;

  Bio::EnsEMBL::Registry->load_all( $regfile );
  $ENS_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $compara, 'compara' );
  $ENS_DBA || ( print( "\n[*DIE] No compara DB set in $regfile\n\n" ) &&
                pod2usage() );
}

open( FH, $SYNTENY_FILE );

my $header = <FH>;
# TODO: Determine the following by parsing the header!
my $REG = "Bio::EnsEMBL::Registry";

#my ($query_species_name, $target_species_name) = ($query_species, $target_species); 
#$query_species_name =~ s/\s+/_/g;
#$target_species_name =~ s/\s+/_/g;

my $meta = {
  'A' => {
    'name'  => "$query_species", #'Zea_mays',
    #'units' => 'band',
  },
  'B' => {
    'name'  => "$target_species", #'Oryza_sativa_japonica',
#    'units' => 'basepair',
  }
};

# Populate some meta data
my $gdba  = $ENS_DBA->get_adaptor('GenomeDB'); # Compara's 'species'
foreach my $key( keys %$meta ){
  my $species = $meta->{$key}->{'name'};
  my $core_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'core') ||
      ( print( "\n[*DIE] No core DB for $species set in registry\n\n" ) &&
        pod2usage() );
  
  $meta->{$key}->{'core_adaptor'}  = $core_adaptor;
  $meta->{$key}->{'slice_adaptor'} = $core_adaptor->get_adaptor('Slice');
  $meta->{$key}->{'genome_db'} = $gdba->fetch_by_registry_name($species) ||
      ( print( "\n[*DIE] Species $species not in compara DB\n\n" )
        && pod2usage() );
}


# Process the synteny file 
my $linenum = 1;
while ( my $line = <FH> ){
  my $syn_meta = parse_line( $line ) || next;
  
  # Get a new synteny_id to use for this line
  my $synteny_id = create_synteny_region($meta->{'A'}->{'genome_db'},
                                         $meta->{'B'}->{'genome_db'});

  foreach my $key( keys %$meta ){
    my $sp_meta     = $meta->{$key};
    my $sp_syn_meta = $syn_meta->{$key};
    
    my $species  = $sp_meta->{name};
    my $seqname  = $sp_syn_meta->{seqname};    
    my $start    = $sp_syn_meta->{start};
    my $end      = $sp_syn_meta->{end};
    my $strand   = $sp_syn_meta->{strand};

    unless ( $strand =~ /^[01]$/){
      $strand = $strand eq '+' ? 1 : 0;
    }

    # Convert all units into basepair
    if( $sp_meta->{units} && $sp_meta->{units} eq 'band' ){
      $start = $start * $BAND_TO_BASEPAIR;
      $end   = $end   * $BAND_TO_BASEPAIR;
    }

    # Convert syntenic region into toplevel
    my $slice = $sp_meta->{slice_adaptor}->fetch_by_region
        ( undef, $seqname, $start, $end, $strand );
    
    $slice || ( print( "\n[*DIE] Sequence $seqname not in core DB for ".
                       "$species\n\n" ) && pod2usage() );
    $slice = ${$slice->project('toplevel')}[0]->to_Slice();
    $seqname = $slice->seq_region_name;
    $start   = $slice->start;
    $end     = $slice->end;
    #$strand  = $slice->strand;
    print "$seqname, $start, $end, $strand\n";

    # Get the DnaFrag corresponding to the top-level slice
    my $dna_frag = get_dna_frag($sp_meta->{genome_db}, $slice );

    # Load the syntenic region into compara
    create_dnafrag_region( $synteny_id, $dna_frag->dbID,
                           $start, $end, $strand );
  }
}


#----------------------------------------------------------------------
# Parses a line from the synteny file, and assigns data to hashref
# Returns undef if problem with line
sub parse_line{
  my $line = shift;
  chomp $line;
  my @parts = split( /\s+/, $line, 8 );
  my @empty = grep{ ! length($_) } @parts[0..7];

  if( @empty ){
    # Sanity check
    print("[WARN] Incomplete line in synteny file: $line\n");
    return undef;
  }
  
  return {
    'A' => { 'seqname' => $parts[0],
             'start'   => $parts[1],
             'end'     => $parts[2],
	     'strand'  => $parts[3], },
    'B' => { 'seqname' => $parts[4],
             'start'   => $parts[5],
             'end'     => $parts[6],
	     'strand'  => $parts[7],},
  };
}

#----------------------------------------------------------------------
# 
sub get_dna_frag{
  my $genome_db = shift || die( "Need a valid GenomeDB object" );
  my $slice     = shift || die( "Need a Slice object" );

  my $seqname = $slice->seq_region_name;

  my $dnafa = $ENS_DBA->get_adaptor('DnaFrag');
  my $dnafrag = $dnafa->fetch_by_GenomeDB_and_name( $genome_db, $seqname );
  $dnafrag && return $dnafrag;

  # Not found in DB - create new
  $dnafrag = Bio::EnsEMBL::Compara::DnaFrag->new
      ( -name              => $seqname,
        -length            => $slice->length,
        -coord_system_name => $slice->coord_system->name, 
        -genome_db         => $genome_db, );
  $dnafa->store( $dnafrag ) if $INSERT;

  return $dnafrag;
}

#----------------------------------------------------------------------
# Inserts a record into the synteny_region table of the ensembl DB and
# returns it's dbID
sub create_synteny_region{
  my @genome_dbs = @_;
  @genome_dbs >= 2 || die( "Need at least two GenomeDBs" );

  my $link_name = "SYNTENY";

  my $mlssa = $ENS_DBA->get_adaptor('MethodLinkSpeciesSet');

  my $mlss = $mlssa->fetch_by_method_link_type_GenomeDBs
      ($link_name,[@genome_dbs]);

  unless( $mlss ){
    $mlss = Bio::EnsEMBL::Compara::MethodLinkSpeciesSet->new
        (
         -method_link_type => $link_name,
         -species_set      => [@genome_dbs],
         );
    $mlssa->store( $mlss );
  }
  
  my $q = qq(
INSERT INTO synteny_region ( method_link_species_set_id )
VALUES (?) );

  my $sth = $ENS_DBA->dbc->prepare( $q );
  $sth->execute( $mlss->dbID ) || die( $sth->errstr );
  return $sth->{'mysql_insertid'};
}

#----------------------------------------------------------------------
#
sub create_dnafrag_region{
  my $synteny_id    = shift || die( "Need a synteny_region_id" );
  my $dnafrag_id    = shift || die( "Need a dnafrag_id" );
  my $dnafrag_start = shift || die( "Need a dnafrag_start" );
  my $dnafrag_end   = shift || die( "Need a dnafrag_end" );
  my $dnafrag_strand   = shift || 0;
  
  my $q = qq(
INSERT INTO dnafrag_region 
       ( synteny_region_id, dnafrag_id, dnafrag_start, dnafrag_end, dnafrag_strand )
VALUES ( ?,?,?,?, ?) );

  my $sth = $ENS_DBA->dbc->prepare( $q );
  $sth->execute( $synteny_id,$dnafrag_id,$dnafrag_start,$dnafrag_end,$dnafrag_strand ) || 
      die( $sth->errstr );
  return $sth->{'mysql_insertid'};
}
