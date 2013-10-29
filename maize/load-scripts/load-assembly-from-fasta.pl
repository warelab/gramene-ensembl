#!/usr/bin/env perl

=head1 NAME

load_assembly_from_fasta.pl - Populates a core Ensembl DB with FASTA data

=head1 SYNOPSIS

perl load_assembly_from_fasta.pl [options] fasta_file(s)

Options:

 -h|--help
 -m|--man
 -e|--ensembl_registry file
 -s|--species species_name
 -a|--assembly_version
 -co|--coord_system name
 -ch|--chunk_size size
 -n|--no_insert

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

B<-a --assembly_version>
  Use this string to indicate the assembly version (def 'unspecified').

B<-co --coord_system>
  The coordinate system appropriate for the sequences (e.g. clone,
  scaffold, contig, chromosome).

B<-ch --chunk_size>
  Split the FASTA entries into chunks of this number of basepairs
  (typically 100000).

B<-n --no_insert>
  Do not make changes to the database. Useful for debug.

=head1 DESCRIPTION

B<This program> 

  Populates the dna, seq_region and coord_system tables of an 'empty'
  core Ensembl DB with data provided by one or more fasta files.

  It also populates various species and assembly-related tables in the
  meta table of the db.

  To create an empty ensembl database, first create an appropriately
  named db from the MySQL client, then load the schema using the file
  in <ENSEMBL_API_ROOT>/ensembl/table.sql

  If the B<--chunk_size> argument is provided, the original sequence
  will be split into seq_region chunks of this size against a
  sequence_level coord_system of 'chunk'. The mapping from the chunks
  to the original B<--coord_system> seq_regions will be loaded into
  the assembly table. This approach is useful for loading whole
  chromosomes where the length of the sequences are larger than the
  size of data that can be loaded into the dna table, or sensibly
  manipulated by the Ensembl API.

  TODO: Add option to chunk sequences where long strings of 'N'
  characters are found.

Maintained by Will Spooner <whs@ebi.ac.uk>

=cut


BEGIN {
    $ENV{'GrameneDir'} ||= '/usr/local/gramene/'; 
    $ENV{'GrameneEnsemblDir'} ||= '/usr/local/gramene_ensembl/'; 
}

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
  # Set the perl libraries. Need a symlink from the gramene-ensembl 
  #   directory to an ensembl API distribution
  $BASEDIR = dirname($Bin);
  unshift @INC, $BASEDIR.'/ensembl-live/ensembl/modules';
  unshift @INC, $BASEDIR.'/ensembl-live/ensembl-compara/modules';
  unshift @INC, $BASEDIR.'/bioperl-HEAD'; # ncomment if BioPerl not in path
}

use Bio::EnsEMBL::Registry; 
use Bio::EnsEMBL::CoordSystem; 
use Bio::EnsEMBL::Slice; 
use Bio::SeqIO::MultiFile;
use Bio::DB::Taxonomy;
use Bio::Species;

use vars qw( $ENS_DBA $SEQ_IO $CS_NAME $CHUNK_SIZE $I $TAXON $VERSION $NOT_SEQ_LEVEL);
BEGIN{  #Argument Processing
  my( $help, $man, $no_insert );
  
  my( $species, $file );
  GetOptions
      (
       "help|?"               => \$help,
       "man"                  => \$man,
       "no_insert"            => \$no_insert,
       "species=s"            => \$species,
       "assembly_version=s"   => \$VERSION,
       "ensembl_registry=s"   => \$file,
       "coord_system=s"       => \$CS_NAME,
       "chunk_size=s"           => \$CHUNK_SIZE,
       "not_seq_level"        => \$NOT_SEQ_LEVEL,
       ) or pod2usage(2);
  pod2usage(-verbose => 2) if $man;
  pod2usage(1) if $help;

  $file    ||= $BASEDIR.'/conf/ensembl.registry';
  $VERSION ||= 'unspecified';
  $CS_NAME || pod2usage( "\nNeed a coord_system name\n" );
  $species || pod2usage( "\nNeed a --species\n" );
  @ARGV    || pod2usage( "\nNeed the path to FASTA file(s)\n" );

  map{
    -e $_ || pod2usage( "\nFile $_ does not exist\n" );
    -r $_ || pod2usage( "\nCannot read $_\n" );
    -f $_ || pod2usage( "\nFile $_ is not plain-text\n" );
    -s $_ || pod2usage( "\nFile $_ is empty\n" );
  } $file, @ARGV;
  
  $I = $no_insert ? 0 : 1;
  unless( $I ){ warn( "[INFO] Evaluation run - no db changes will be made\n")}

  # Load the ensembl file
  Bio::EnsEMBL::Registry->load_all( $file );
  $ENS_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'core' );
  $ENS_DBA || pod2usage( "\nNo core DB for $species set in $file\n" );

  # Load the Fasta file
  warn( "[INFO] FASTA files: " . join( "\n".(' 'x18), @ARGV ) . "\n" );
  $SEQ_IO = new Bio::SeqIO::MultiFile( -format  => "fasta",
                                       -files   => [@ARGV]);

  # Get taxonomy info for the species
  my $species_binomial = Bio::EnsEMBL::Registry->get_alias($species);
  $species_binomial =~ s/_/ /g;
#  $species_binomial =~ s/scaffold//g;
  my $db = Bio::DB::Taxonomy->new(-source => 'entrez' );
  my $taxonid = $db->get_taxonid( $species_binomial )
      || pod2usage( "\nNo NCBI taxon id for $species_binomial".
                    "\nTry using a more recent bioperl version\n" );
  $TAXON = $db->get_Taxonomy_Node( -name => $species_binomial )
      || pod2usage( "\nNo NCBI taxonomy for $species_binomial\n" );
}

#==================== MAIN ====================

# Deal with the taxon and assembly version
# &store_taxon_meta( $TAXON );
my $mc = $ENS_DBA->get_MetaContainer;
warn( "[INFO] adding assembly.default $VERSION to meta table\n" );
if( $I ){ $mc->store_key_value('assembly.default', $VERSION) }

# Get the coord_systems
my( $seq_cs, $chunk_cs ) = &get_coord_systems();
my $slice_adaptor = $ENS_DBA->get_adaptor('Slice');

while( my $seq = $SEQ_IO->next_seq ){

  my $seq_id     = $seq->id; # TODO - use filename if ID is missing

  if($CS_NAME =~ /chromosome/i){
    next unless ($seq_id =~ s/chr(omosome)?_?0*//i);# Strip 'chromsome' prefix from ID
    
  }else{
    next if ($seq_id =~ s/chr(omosome)?_?0*//i);
  }
  #$seq_id =~ s/chr(omosome)?_?0*//i; # Strip 'chromsome' prefix from ID
  $seq->id( $seq_id );

  my $seq_length = $seq->length;

  warn( "[INFO] Processing $CS_NAME $seq_id ($seq_length bp)\n" );
  my $seq_slice = &new_slice( $seq_cs, $seq );
  if( $CHUNK_SIZE ){
    for( my $start=1; $start<=$seq_length; $start+=$CHUNK_SIZE ){
      my $end = $start + $CHUNK_SIZE - 1;
      if( $end > $seq_length ){ $end = $seq_length }
      my $chunk_id = "${seq_id}_${start}..${end}";
      warn( "[INFO]   Processing chunk $chunk_id\n" );
      my $chunk_seq = Bio::PrimarySeq->new
          (
           -id => $chunk_id,
           -seq => $seq->subseq( $start, $end ),
           );
      my $chunk_slice = &new_slice( $chunk_cs, $chunk_seq ); 
      my $asm_slice = $seq_slice->sub_Slice( $start, $end );
      if( $I ){
        $slice_adaptor->store_assembly( $asm_slice, $chunk_slice );
      }
    }
  }
}


#======================================================================
# Takes a Bio::Taxonomy::Node and populates the species
# entries in the meta table
sub store_taxon_meta{
  my $taxon = shift || die( "Need a Bio::Taxonomy::Node" );
  my $mc = $ENS_DBA->get_MetaContainer;

  my $taxid = $taxon->object_id;
  warn( "[INFO] adding species.taxonomy_id $taxid to meta\n" );
  if( $I ){ $mc->store_key_value('species.taxonomy_id', $taxid) }

  if( my $name = $taxon->common_name ){
    warn( "[INFO] adding species.common_name $name to meta\n" );
    if( $I ){ $mc->store_key_value('species.common_name', $name) }
  }

  foreach my $name( $taxon->classification ){
    my $key = 'species.classification';
    warn( "[INFO] adding $key $name to meta\n" );
    if( $I ){ $mc->store_key_value($key, $name) }
  }
  return;
}

#----------------------------------------------------------------------
# Creates the appropriate coord systems in the DB
sub get_coord_systems{

  my $csa = $ENS_DBA->get_adaptor('CoordSystem');

  my $seq_cs ;
  my $seq_cs_exist;

  $seq_cs = $csa->fetch_by_name($CS_NAME);

  if($seq_cs){

    $seq_cs_exist = 1;
  }else{
  
    $seq_cs = Bio::EnsEMBL::CoordSystem->new
      (
       -name    => $CS_NAME,
       -rank    => $CS_NAME =~ /chromosome/i ? 1 : 2,
       -default => 1,
       -version => $VERSION, 
       $CHUNK_SIZE ? () : 
       ( $NOT_SEQ_LEVEL ? () : -sequence_level=>1 ),
       );
  }

  my $chunk_cs;
  my $chunk_cs_exist;
  if( $CHUNK_SIZE ){
    $chunk_cs = $csa->fetch_by_name('chunk');
    
    if($chunk_cs){
      $chunk_cs_exist = 1;
    }else{
      $chunk_cs = Bio::EnsEMBL::CoordSystem->new
	(
	 -name           => 'chunk',
	 -rank           => 3,
	 -default        => 1,
	 -version => $VERSION,
	 -sequence_level => 1,
	);
    }
  } 
  
  warn( "[INFO] Storing coord_system: $CS_NAME\n" );
  if( $I ){ $csa->store($seq_cs) unless $seq_cs_exist; }
  if( $chunk_cs ){
    print( "[INFO] Storing coord_system: chunk ($CHUNK_SIZE bp)\n" );
    if( $I ){
      $csa->store($chunk_cs) unless $chunk_cs_exist;
      $csa->store_mapping_path($seq_cs, $chunk_cs );
    }
  }
  return( $seq_cs, $chunk_cs );
}

#----------------------------------------------------------------------
# Creates a new slice based on the given cs_name, name, start and end
our $SLA;
sub new_slice{
  my $cs = shift || die( "Need a coord system" );
  my $bioseq = shift || die( "Need a Bio::Seq" );

#printf "seq_region_name=%s, len=%d\n", $bioseq->id, $bioseq->length;

  my $slice = Bio::EnsEMBL::Slice->new
      (
       -coord_system      => $cs,
       -seq_region_name   => $bioseq->id,
       -start             => 1,
       -end               => $bioseq->length,
       -strand            => 1,
       -version           => 1,
       );

  my $seqref;
  if( $cs->is_sequence_level ){ $seqref = \($bioseq->seq) }

  if( $I ){
    $SLA ||= $ENS_DBA->get_adaptor('Slice');
    $SLA->store( $slice, ( $seqref ? $seqref : () ) )
  }
  return $slice;
}

#======================================================================
1;
__END__

