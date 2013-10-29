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
 -r|--rank
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

B<-r --rank>
  The rank of this coordinate system 

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
    $ENV{'GrameneEnsemblDir'} ||= '/usr/local/ensembl-live/'; 
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

use vars qw( $ENS_DBA $SEQ_IO $CS_NAME $CHUNK_SIZE $I $TAXON $VERSION 
$NOT_SEQ_LEVEL $RANK $FORCE_SEQ);
BEGIN{  #Argument Processing
  my( $help, $man, $no_insert );
  
  my( $species, $file );
  GetOptions
      (
       "help|?"               => \$help,
       "man"                  => \$man,

       "species=s"            => \$species,

       ) or pod2usage(2);
  pod2usage(-verbose => 2) if $man;
  pod2usage(1) if $help;

  $species || pod2usage( "\nNeed a --species\n" );


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



#======================================================================
# Takes a Bio::Taxonomy::Node and populates the species
# entries in the meta table
sub store_taxon_meta{
  my $taxon = shift || die( "Need a Bio::Taxonomy::Node" );


  my $taxid = $taxon->object_id;
  warn( "[INFO] adding species.taxonomy_id $taxid to meta\n" );

  if( my $name = $taxon->common_name ){
    warn( "[INFO] adding species.common_name $name to meta\n" );
  }

  use Bio::Tree::Tree;
  my $tree_functions = Bio::Tree::Tree->new();
  

  #use Data::Dumper Dumper;
  foreach my $node( reverse $tree_functions->get_lineage_nodes($taxon) ){
    #warn( Dumper ($node));
    my $name = $node->node_name;
    my $key = 'species.classification';
    warn( "[INFO] adding $key $name to meta\n" );

  }
  return;
}



#======================================================================
1;
__END__

