#!/bin/env/perl -w

=pod

=head1 NAME

perl dump_tree_homologs_by_species.pl - Dumps homologs/paralogs for each gene

=head1 SYNOPSIS

  perl dump_tree_id.pl [options] > resultfile.dat

Options:

  -h|--help             Show brief help and exit.
  -m|--man              Show detailed help
  -e|--ensembl_registry Path to Ensembl registry file.
  -s|--species          Species key in Ensembl registry file.
  -t|--taxonomy_id      Dump related genes from these species tax ids
  -sl|--slice_name      Name of a slice to dump from NOT IMPLEMENTED

=head1 OPTIONS

B<-h|--help>
  Print a brief help message and exits.

B<-m|--man>
  Print man page and exit

B<-e|-ensembl_registry>
  Use this Ensembl registry file for database connection info.
  Default is <ENSEMBLHOME>/conf/ensembl.registry

B<-s|--species>
  Reference genes from this species in the registry file [REQUIRED].

B<-t|--taxonomy_id>
  Dump relationships between reference genes and orthologs/paralogs in
  these related species, identified by NCBI taxonomy ID (can be more than one).
  Defaults to all species.
    
B<-sl|--slice_name>
  NOT YET IMPLEMENTED
  Dump reference genes only from the slice identified by this name. 
  See the Ensembl documentation for the Bio::EnsEMBL::Slice object
  for info on the format. Defaults to whole genome.

=head1 DESCRIPTION

Script that iterates through each gene in the given species and dumps
information about their orthologs/paralogs. 

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

use vars qw( $BASEDIR );
BEGIN{
  # Set the perl libraries 
  $BASEDIR = dirname($Bin);
  -d $BASEDIR.'/ensembl-live'
    || die( "\n[*DIE] Need $BASEDIR/ensembl-live symlinked to ensembl\n" );
  unshift @INC, $BASEDIR.'/ensembl-live/ensembl/modules';
  unshift @INC, $BASEDIR.'/ensembl-live/ensembl-compara/modules';
}

use Bio::EnsEMBL::Registry;

our $ENS_DBA;
our $CMP_DBA;
our @REL_SPECIES;
our $SLICE;
BEGIN{  #Argument Processing
  my $help=0;
  my $man=0;
  my( $species, $reg, $slice_name );
  GetOptions
      ( 
        "help|?"             => \$help,
        "man"                => \$man,
        "s|species=s"        => \$species,
        "ensembl_registry=s" => \$reg,
        "taxonomy_id=s"      => \@REL_SPECIES,
        "slice_name=s"       => \$slice_name,
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
      || pod2usage("\n[*DIE] No core DB for $species set in $reg\n" );
  $CMP_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( 'compara','compara')
      || pod2usage("\n[*DIE] No compara DB set in $reg\n" );
  my $sla = $ENS_DBA->get_adaptor('Slice')
      || pod2usage( "[*DIE] Cannot ENS_DBA->get_adaptor('Slice')" );
  if( $slice_name ){
    $SLICE = $sla->fetch_by_name($slice_name)
        || pod2usage( "[*DIE] Cannot SLICE->fetch_by_name($slice_name)\n" );
    pod2usage( "[*DIE] Sorry, slice arg not yet implemented\n" );
  }
}


MAIN:{
 
  # print the header
  print join( "\t",
              'gene_stable_id',
#              'sequence_type',
#              'sequence_name',
#              'gene_bp_start_in_sequence',
#              'gene_bp_end_in_sequence',
              'root_node_id',
              'root_taxon_name',
              'speciation_node_id',
              'speciation_taxon_name',
              'speciation_mode',
              'homol_gene_stable_id',
              'homol_taxon_name',
              'homol_taxon_id',
              ) . "\n";

  # Create the adaptors that we will be using to fetch data from the database
  my $gene_adaptor = $ENS_DBA->get_adaptor('Gene')
      || die( "[*DIE] Cannot ENS_DBA->get_adaptor('Gene')" );
  my $member_adaptor = $CMP_DBA->get_adaptor('Member')
      || die( "[*DIE] Cannot CMP_DBA->get_adaptor('Member')" );
  my $tree_adaptor  = $CMP_DBA->get_adaptor('ProteinTree');


  # These are the taxon IDs of homologous species to include
  my %search_taxon_ids = map{ $_ => 1 } @REL_SPECIES; 

  # Some counters
  my( $num_genes, $num_coding, $num_tree_genes, %tree_nodes ) = (0,0,0,());

  # Degug
  my $db = $ENS_DBA->dbc->dbname;
  warn( "[INFO] Collecting all genes from $db \n" );

  # Get a list of all gene IDs, and loop through them
  my $gene_id_arrayref = $gene_adaptor->list_stable_ids;
  my $total_genes = scalar @$gene_id_arrayref;

  warn( "[INFO] Collected genes $total_genes; processing...\n" );

  # loop through them
  while( my $id = shift @$gene_id_arrayref ){
    $num_genes ++; #COUNT

    # Get the gene object, and skip unless protein coding
    my $gene = $gene_adaptor->fetch_by_stable_id( $id )
        || die( "[*DIE] Cannot fetch_by_stable_id gene $id" );
    next unless $gene->biotype eq 'protein_coding';
    $num_coding ++; # COUNT

    # Get the subroot for the gene tree
    my $member  
        =  $member_adaptor->fetch_by_source_stable_id('ENSEMBLGENE',$id)
        || ( printf("%s\n", $id ) && next );

    my $gene_tree
        = $tree_adaptor->fetch_by_Member_root_id($member, 0)
        || ( printf("%s\n", $id ) && next );

    my $root_node   = $gene_tree->subroot;
    my $root_node_id = $root_node->node_id;
    my $root_taxon   = $root_node->get_tagvalue('taxon_name') || '';
    $num_tree_genes                ++; #COUNT
    $tree_nodes{$root_node_id} ++; #COUNT
    
    my @gene_data = ( $id,
                      $root_node_id,
                      $root_taxon );
    
    my $gene_node;
    my @homol_nodes;
    foreach my $leaf_node( @{$root_node->get_all_leaves} ){
      if( $leaf_node->gene_member->stable_id eq $id ){ #this gene
        $gene_node = $leaf_node;
        next;
      }
      elsif( %search_taxon_ids and 
          ! $search_taxon_ids{$leaf_node->genome_db->taxon_id} ){ next }
      else{
        push @homol_nodes, $leaf_node;
      }
    }
    
    foreach my $homol_node( @homol_nodes ){
      my $common_node = $gene_node->find_first_shared_ancestor($homol_node);
      my $speciation_node_id = $common_node->node_id;
      my $speciation_taxon   = $common_node->get_tagvalue('taxon_name') || '';
      my $speciation_mode = $common_node->get_tagvalue('Duplication') 
          ? ( $common_node->get_tagvalue('dubious_duplication') 
              ? 'Dubious_duplication' : 'Duplication' )
          : 'Speciation';
      my $homol_gene_id = $homol_node->gene_member->stable_id;
      my $homol_taxon = $homol_node->genome_db->name || '';
      my $homol_taxon_id = $homol_node->genome_db->taxon_id;
      # Print the data
      print join( "\t", 
                  @gene_data,
                  $speciation_node_id,
                  $speciation_taxon,
                  $speciation_mode,
                  $homol_gene_id,
                  $homol_taxon,
                  $homol_taxon_id,
                  ) . "\n";
      

    }

    unless( @homol_nodes ){
      # Print the data
      print join( "\t", @gene_data) . "\n";
    }

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
