#!/usr/bin/env perl

=head1 NAME

load_assembly_from_gamexml.pl - Populates a core Ensembl DB with gamexml data

=head1 SYNOPSIS

perl load_assembly_from_gamexml.pl [options] gamexml_file(s)

Options:

 -h|--help
 -m|--man
 -e|--ensembl_registry file
 -s|--species species_name
 -a|--assembly_version
 -co|--coord_system name
 
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


B<-n --no_insert>
  Do not make changes to the database. Useful for debug.

=head1 DESCRIPTION

B<This program> 

  Populates the dna, seq_region and coord_system tables of an 'empty'
  core Ensembl DB with data provided by one or more gamexml files.

  It also populates various species and assembly-related tables in the
  meta table of the db.

  To create an empty ensembl database, first create an appropriately
  named db from the MySQL client, then load the schema using the file
  in <ENSEMBL_API_ROOT>/ensembl/table.sql

  The seq_level would be the chunk, corresponding to the gamexml segments
  The upper level coord_system would be defined from command line arg. The
  possible ones could be chromosome, superscaffold.


  TODO: Add option to chunk sequences where long strings of 'N'
  characters are found.

Modified from load_assembly_from_fasta.pl by Sharon Wei (weix@cshl.edu)

=cut


BEGIN {
    $ENV{'GrameneDir'} = '/usr/local/gramene/'; 
    $ENV{'GrameneEnsemblDir'} = '/usr/local/ensembl-live/'; 
}


use lib map {"$ENV{'GrameneEnsemblDir'}/$_"} 
  qw( ensembl/modules
      ensembl-compara/modules
      );

use lib '/usr/local/lib/perl5/site_perl/5.8.8';

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper qw(Dumper); # For debug 
use DBI;
use FindBin qw( $Bin );
use File::Basename qw( dirname basename );
use Bio::SeqIO;

use Bio::EnsEMBL::Registry; 
use Bio::EnsEMBL::CoordSystem; 
use Bio::EnsEMBL::Slice; 
use Bio::SeqIO::MultiFile;
use Bio::DB::Taxonomy;
use Bio::Species;

use vars qw( $ENS_DBA $CS_NAME $CHUNK_SIZE $I $TAXON $VERSION 
$NOT_SEQ_LEVEL $RANK $FORCE_SEQ);
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
       "rank=i"               => \$RANK,
       "chunk_size=s"           => \$CHUNK_SIZE,
       "not_seq_level"        => \$NOT_SEQ_LEVEL,
       "force_seq"            => \$FORCE_SEQ,
       ) or pod2usage(2);
  pod2usage(-verbose => 2) if $man;
  pod2usage(1) if $help;

  $file    or pod2usage( "\nNeed a registry file≈\n" );;
  $VERSION ||= 'unspecified';
  $CS_NAME || pod2usage( "\nNeed a coord_system name\n" );
  $species || pod2usage( "\nNeed a --species\n" );
  @ARGV    || pod2usage( "\nNeed the path to GAMEXML file(s)\n" );

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
  warn( "[INFO] GAMEXML files: " . join( "\n".(' 'x18), @ARGV ) . "\n" );
 


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
&store_taxon_meta( $TAXON );

my $mc = $ENS_DBA->get_MetaContainer;
warn( "[INFO] adding assembly.default $VERSION to meta table\n" );
if( $I ){ $mc->store_key_value('assembly.default', $VERSION) }

# Get the coord_systems
my( $asm_cs, $chunk_cs ) = &get_coord_systems();
my $slice_adaptor = $ENS_DBA->get_adaptor('Slice');

#get chr info
my %chr;

for my $file( @ARGV ){
  
  my $SEQ_IO = Bio::SeqIO->new('-format' => 'game',
			       '-file'   => $file );

  my $file_basename = basename ($file, ".xml" );
  
  my ($asm_seq_region_name, $asm_start, $asm_end) = 
    $file_basename =~ /_chrom0*(\d+)_(\d+)_(\d+)$/i;
  
  my $seq = $SEQ_IO->next_seq ;
    
  my $seq_id  = $seq->id; 
    
  if( $seq_id !~ /chr(omosome)?_?0*(\d+)_?/i ){
    die ("Unexpected seq_id $seq_id");
    #next;
  }
   
  #                     genomic_start_coord, genomic_end_coord, seq_string
  push @{$chr{$asm_seq_region_name}}, [$asm_start, $asm_end, $seq->seq, $file_basename];
  
}

for my $chr(sort keys %chr){
    
  #print "chr=$chr, $chr{$chr}-

  $chr{$chr} = [sort {$a->[0] <=> $b->[0]} @{$chr{$chr}}];
}

exit unless $I;


  
for my $chr(sort keys %chr){
  
  my ($chr_start, $chr_end) = ($chr{$chr}->[0]->[0], $chr{$chr}->[-1]->[1]);
  my $chunk_cnt = scalar @{$chr{$chr}};
  
  print "chr=$chr, $chr_start, $chr_end\n";
  
  die "chr start ($chr_start) is not 1" if $chr_start != 1;
  
  
  unless( $slice_adaptor->fetch_by_region( $CS_NAME, $chr ) ){
    
    
    my $chr_slice = Bio::EnsEMBL::Slice->new
      (
       -coord_system      => $asm_cs,
       -seq_region_name   => $chr,
       -start             => 1,
       -end               => $chr_end,
       -strand            => 1,
       -version           => 1,
      );
    
    $slice_adaptor->store( $chr_slice );
  }
  
  my $pre_start = 0;
  my $pre_end   = 0;

  for (my $i=0; $i<$chunk_cnt; $i++){
      
    my ( $cur_start, $cur_end, $seq_str, $file_basename ) = @{$chr{$chr}->[$i]};
    my $seq_length = length( $seq_str);
    if( $seq_length != ($cur_end-$cur_start+1)){
      die ("unmatch seq length, $seq_length != ($cur_end-$cur_start+1)\n");
      next;
    }

    warn( "[INFO] Processing $CS_NAME  ($seq_length bp) - $file_basename\n" );

    unless( $slice_adaptor->fetch_by_region( 'chunk', $file_basename ) ){
      my $chunk_slice = Bio::EnsEMBL::Slice->new
	(
	 -coord_system      => $chunk_cs,
	 -seq_region_name   => $file_basename,
	 -start             => 1,
	 -end               => $seq_length,
	 -strand            => 1,
	 -version           => 1,
	);
      $slice_adaptor->store( $chunk_slice, \$seq_str );
    }
    
    my ($asm_slice,  $asm_start, $asm_end);
    my ($comp_slice, $comp_start, $comp_end);
    
    if( $i > 0 && $cur_start <= $pre_end ){
	my $overlap   = $pre_end - $cur_start + 1;
	
	$asm_start = $cur_start + $overlap;
	$asm_end   = $cur_end;
	
	$comp_start = 1 + $overlap;
	$comp_end   = $seq_length;
	
      }else{
	
	($asm_start, $asm_end) = ( $cur_start, $cur_end);
      
	($comp_start, $comp_end) = ( 1, $seq_length);
      
      }
    
    $asm_slice = $slice_adaptor->fetch_by_region( 
						 $CS_NAME, 
						 $chr,
						 $asm_start,
						 $asm_end,
						 1);
    $comp_slice = $slice_adaptor->fetch_by_region( 
						  'chunk', 
						  $file_basename,
						  $comp_start,
						  $comp_end,
						  1);
    
    $slice_adaptor->store_assembly( $asm_slice, $comp_slice );
    
    ($pre_start, $pre_end) = ($cur_start, $cur_end);
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

  use Bio::Tree::Tree;
  my $tree_functions = Bio::Tree::Tree->new();
  

  #use Data::Dumper Dumper;
  foreach my $node( $tree_functions->get_lineage_nodes($taxon) ){
    #warn( Dumper ($node));
    my $name = $node->node_name;
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

  my $asm_cs = $csa->fetch_by_name($CS_NAME);
  
  my $asm_cs_store;
  
  unless ($asm_cs){
    $asm_cs = Bio::EnsEMBL::CoordSystem->new
      (
       -name    => $CS_NAME,
       -rank    => 1
       -default => 1,
       -version => $VERSION, 
       );
    $asm_cs_store = 1;
    
  }

  my $chunk_cs;
  my $chunk_cs_store;
  
  $chunk_cs = $csa->fetch_by_name('chunk');
    
  unless($chunk_cs){
   
    $chunk_cs = Bio::EnsEMBL::CoordSystem->new
      (
       -name           => 'chunk',
       -rank           => 2,
       -default        => 1,
       -version        => $VERSION,
       -sequence_level => 1,
      );
    $chunk_cs_store = 1;
  
  } 
  
  warn( "[INFO] Storing coord_system: $CS_NAME\n" );
  if( $I ){ 

    $csa->store($asm_cs) if $asm_cs_store;
    $csa->store($chunk_cs) if $chunk_cs_store;
    $csa->store_mapping_path($asm_cs, $chunk_cs );
  }
	
  return( $asm_cs, $chunk_cs );
}

#----------------------------------------------------------------------
# Creates a new slice based on the given cs_name, name, start and end
our $SLA;
sub new_slice{
  my $cs = shift || die( "Need a coord system" );
  my $bioseq = shift || die( "Need a Bio::Seq" );

#printf "seq_region_name=%s, len=%d\n", $bioseq->id, $bioseq->length;

  #gi|84095186|dbj|AP008219.1|
  my $seq_region_name = $bioseq->display_id;
  if($seq_region_name =~ /gi\|\d+\|\w+\|(\w+)/i){
    $seq_region_name = $1;
  }

  my $slice = Bio::EnsEMBL::Slice->new
      (
       -coord_system      => $cs,
       -seq_region_name   => $seq_region_name,
       -start             => 1,
       -end               => $bioseq->length,
       -strand            => 1,
       -version           => 1,
       );

  my $seqref;
  if( $cs->is_sequence_level || $FORCE_SEQ ){ $seqref = \($bioseq->seq) }

  if( $I ){
    $SLA ||= $ENS_DBA->get_adaptor('Slice');
    $SLA->store( $slice, ( $seqref ? $seqref : () ) )
  }
  return $slice;
}

#======================================================================
1;
__END__

