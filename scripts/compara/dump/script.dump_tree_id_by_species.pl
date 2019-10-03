#!/usr/bin/env perl

=pod

=head1 NAME

perl dump_tree_id_by_species.pl - Dumps the protein tree ID for each gene

=head1 SYNOPSIS

  perl dump_tree_id.pl [options] > resultfile.dat

Options:

  -h|--help             Show brief help and exit.
  -m|--man              Show detailed help
  -e|--ensembl_registry Path to Ensembl registry file.
  -s|--species          Species key in Ensembl registry file.

=head1 OPTIONS

B<-h|--help>
  Print a brief help message and exits.

B<-m|--man>
  Print man page and exit

B<-e|-ensembl_registry>
  Use this Ensembl registry file for database connection info.
  Default is <ENSEMBLHOME>/conf/ensembl.registry

B<-s|--species>
  Use this species entry from the registry file [REQUIRED].

=head1 DESCRIPTION

Script that iterates through each gene in the given species and dumps
the subroot ID of any protein tree alongside the gene genomic
coordinates

B<The Ensembl Registry>

  The database connection details for both Ensembl species and compara
  databases must be configured in an ensembl registry file.

  An example registry file would be;
  ---
  use Bio::EnsEMBL::DBSQL::DBAdaptor; 
  use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
  use Bio::EnsEMBL::Registry;
  # The Ensembl core database
  Bio::EnsEMBL::DBSQL::DBAdaptor->new
  ( '-species' => 'Oryza_sativa', 
    '-group'   => 'core', 
    '-dbname'  => 'Oryza_sativa_japonica_core_48_28', );
  Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new
  ( '-species' => 'compara',
    '-group'   => 'compara',
    '-dbname'  => 'ensembl_compara_48_28', );
  ---

TODO: Complete this section
                   
Maintained by Will Spooner <whs@ebi.ac.uk>

=cut


use strict;
use warnings;
use Data::Dumper qw(Dumper);
use File::Basename;
use FindBin qw( $Bin );
use Pod::Usage;
use Getopt::Long;
use IO::File;

use vars qw( $BASEDIR $INTERNAL );
BEGIN{
  # Set the perl libraries 
    $BASEDIR = $ENV{COMPARA_BASEDIR} || $ENV{ENS_ROOT}; #dirname( dirname( dirname($Bin) ) );
  -d $BASEDIR
      || die( "\n[*DIE] Need $BASEDIR with ensembl modules\n" );
    unshift @INC, $BASEDIR.'/ensembl/modules';
    unshift @INC, $BASEDIR.'/ensembl-compara/modules';
}

use Bio::EnsEMBL::Registry;

our $ENS_DBA;
our $CMP_DBA;
BEGIN{  #Argument Processing
  my $help=0;
  my $man=0;
  my( $species, $old_species, $reg, $logic_name, $no_insert );
  GetOptions
      ( 
        "help|?"             => \$help,
        "man"                => \$man,
        "species=s"          => \$species,
        "ensembl_registry=s" => \$reg,
        )
      or pod2usage(2);
  pod2usage(-verbose => 2) if $man;
  pod2usage(1) if $help;

  $reg        ||= $BASEDIR.'/conf/ensembl.registry';
  map{
    -e $_ || pod2usage( "\n[*DIE]File $_ does not exist\n" );
    -r $_ || pod2usage( "\n[*DIE]Cannot read $_\n" );
    -f $_ || pod2usage( "\n[*DIE]File $_ is not plain-text\n" );
    -s $_ || pod2usage( "\n[*DIE]File $_ is empty\n" );
  } $reg;
  
  # Load the ensembl file
  Bio::EnsEMBL::Registry->load_all( $reg );
  unless( $species ){
    $species    || pod2usage("\n[*DIE] Need a --species\n");
  }
  $ENS_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'core' )
      || pod2usage("\n[*DIE]No core DB for $species set in $reg\n" );
  $CMP_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( 'compara','compara')
      || pod2usage("\n[*DIE]No compara DB set in $reg\n" );

  print "#species:$species\n";
}


MAIN:{

  # print the header
  print join( "\t",
              '#gid',
#              'sequence_type',
              'ref_name',
              'start',
              'end',
	      'strand',
              'root_node_id',
	      'tree_stable_id',
              'root_taxon' ) . "\n";

  # Create the adaptors that we will be using to fetch data from the database
  my $gene_adaptor = $ENS_DBA->get_adaptor('Gene')
      || die( "[*DIE] Cannot ENS_DBA->get_adaptor('Gene')" );
  my $member_adaptor = $CMP_DBA->get_adaptor('GeneMember')
      || die( "[*DIE] Cannot CMP_DBA->get_adaptor('GeneMember')" );
  my $tree_adaptor  = $CMP_DBA->get_adaptor('GeneTree');

  # Some counters
  my( $num_genes, $num_coding, $num_tree_genes, %tree_nodes ) = (0,0,0,());

  # Degug
  my $db = $ENS_DBA->dbc->dbname;
  warn( "[INFO] Collecting all protein coding genes from $db \n" );

  # Get a list of all gene IDs, and loop through them
  my $gene_id_arrayref = $gene_adaptor->list_stable_ids;
  my $total_genes = scalar @$gene_id_arrayref;

  warn( "[INFO] Collected genes $total_genes; processing...\n" );

  # loop through them
  while( my $id = shift @$gene_id_arrayref ){
    $num_genes ++; #COUNT

    # Get the gene object, and skip unless protein coding
    my $gene;

    eval { 
	$gene = $gene_adaptor->fetch_by_stable_id( $id )
	};
     if($@){
	warn( "[*ERROR] Cannot fetch_by_stable_id gene $id, $@" );
	next;
	};
    next unless $gene->biotype eq 'protein_coding';
    $num_coding ++; # COUNT

    # Get the subnode for the gene tree
    my $tree_id;
    my $gene_member  
        =  $member_adaptor->fetch_by_stable_id($id);  #fetch_by_source_stable_id('ENSEMBLGENE',$id);
        #|| die( "Cannot get ENSEMBLGENE compara member for $id" );

    next unless $gene_member;
    
    my $gene_tree;
    if( $gene_tree 
        = $tree_adaptor->fetch_default_for_Member($gene_member) ){
      $tree_id   = $gene_tree->root->node_id;
      my $tree_stable_id = $gene_tree->stable_id;
      
      $num_tree_genes                ++; #COUNT
      $tree_nodes{$tree_id} ++; #COUNT
    }

    # Make sure the gene is on the longest assembled sequence
    $gene->project('toplevel');

    # Print the data
    print join( "\t",
                $id,
#                $gene->coord_system_name,
                $gene->seq_region_name,
                $gene->seq_region_start,
                $gene->seq_region_end,
		$gene->seq_region_strand,
                $gene_tree
                ? ( $gene_tree->root->node_id,	    
		    $gene_tree->stable_id, #does not exist for oge trees
                    $gene_tree->root->taxonomy_level() || '' )
                : ( '','' ) );
    print "\n";

    $gene_tree->release_tree() if $gene_tree;
    # The following lists all node annotations if you want to add extra
    # my @tags = $subroot->get_all_tags;
    # warn join( ", ", @tags );

    # Update progress each 1000 genes
    if( $num_genes/1000-int($num_genes/1000) == 0 ){
      warn( "[INFO] processed $num_genes of $total_genes\n" );
    }
    
  }

  # Debug
  warn sprintf( "[INFO] processed %d coding genes; %d in %d unique trees\n",
                $num_coding, $num_tree_genes, scalar( keys( %tree_nodes ) ) ); 

}

exit;
