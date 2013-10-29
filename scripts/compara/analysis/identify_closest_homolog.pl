#!/bin/env perl -w

# For each gene;
# - Gene stable_id
# - Gene location
# - Best match homolog (use distance to homolog if tie break)
# - Homolog location
# - distance to homolog
# - Alignment length bp (use CDS)
# - Alignment missmatches with homolog
# - Alignment overhang bp
# - Tree ID
# - Tree size
# - Root species
# - Maize-specific tree size
# - All genes with same alignment score.

# This script has to be run from brie to get the paths correct

use strict;
use Data::Dumper qw(Dumper);

# Some constants
our $REGFILE = '/usr/local/gramene/conf/ensembl.registry';
our $SPECIES = 'Zea_mays';
our $COMPARA = 'compara';
our $METHOD_LINK = 'ENSEMBL_PARALOGUES';
our $MAX_DISTANCE= 1000000;
our $OUTFILE = 'Zea_mays_gene_closest_homolog.tsv';
our @FIELDS  = qw( gene_stable_id
                   gene_location
                   homolog_stable_id
                   homolog_location
                   distance_to_homolog_bp
                   gene_length_aa
                   homolog_length_aa
                   alignment_length_aa
                   alignment_missmatches
                   alignment_overhang
                   tree_node_id
                   tree_taxon_name
                   tree_leaf_count
                   species_specific_leaf_count
                   all_homologs
                   );

BEGIN{

  our $BASEDIR = '/usr/local/ensembl-live';

  # Set the perl libraries 
  -d $BASEDIR
    || die( "\n[*DIE] Need a path to the ensembl distribution\n" );
  unshift @INC, $BASEDIR.'/ensembl-live/ensembl/modules';
  unshift @INC, $BASEDIR.'/ensembl-live/ensembl-compara/modules';
}
use Bio::EnsEMBL::Registry;

unless( -f $REGFILE ){
  die( "\n[*DIE] Cannot find ensembl registry $REGFILE\n" );
}

Bio::EnsEMBL::Registry->load_all( $REGFILE );
our $ENS_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $SPECIES, 'core' )
    || die("\n[*DIE] No core DB for $SPECIES set in $REGFILE\n" );
our $CMP_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $COMPARA,'compara')
    || die("\n[*DIE] No compara DB for $COMPARA set in $REGFILE\n" );

my $meta_container = $ENS_DBA->get_adaptor('MetaContainer');
my $gene_adaptor   = $ENS_DBA->get_adaptor('Gene');
my $slice_adaptor  = $ENS_DBA->get_adaptor('Slice');

my $genomedb_adaptor = $CMP_DBA->get_adaptor('GenomeDB');
my $mlss_adaptor     = $CMP_DBA->get_adaptor('MethodLinkSpeciesSet');
my $homology_adaptor = $CMP_DBA->get_adaptor('Homology');
my $tree_adaptor     = $CMP_DBA->get_adaptor('ProteinTree');
our $MEMBER_ADAPTOR   = $CMP_DBA->get_adaptor('Member');

my ($tax_id) = @{$meta_container->list_value_by_key('species.taxonomy_id')};
$tax_id || die( "\n[*DIE] No meta species.taxonomy_id in core DB!\n" );
my( $genome_db ) = ( grep{ $_->taxon_id == $tax_id } 
                     @{$genomedb_adaptor->fetch_all} );
$genome_db || die( "\n[*DIE] No taxon_id $tax_id in compara.\n" );

#my $mlss = $mlss_adaptor->fetch_by_method_link_type_GenomeDBs
#    ( $METHOD_LINK, [$genome_db] )
#    || die( "\n[*DIE] No mlss for $METHOD_LINK and taxon $tax_id.\n" );

# Open the output file
open( OUT, "> $OUTFILE" ) 
    || die("\n[*DIE] Could not open $OUTFILE for write $!\n");
print OUT join( "\t", @FIELDS ). "\n";

#----------------------------------------------------------------------
# This is where the meat of the process is carried out
my %seen_genes;
my %identical_genes;
my %overhanging_genes;
my %single_missmatch_genes;
my %untransformable;
my $dup_count = 0;
MAIN:{
  
  my %genes; # Keep a list of prcessed genes
  # Loop for each top-level slice rather than each gene to keep mem down
  my $slices_ref = $slice_adaptor->fetch_all('toplevel');
  while( my $slice = shift @$slices_ref ){

  #my $genes_ref = $gene_adaptor->fetch_all();
  #my $genes_ref = [$gene_adaptor->fetch_by_stable_id('GRMZM2G081840')];
  my $genes_ref = $slice->get_all_Genes();

  while( my $gene = shift @$genes_ref  ){
    print STDERR '.';
    my %record;
    #$gene = $gene->transform('toplevel');
    my $gene_stable_id = $gene->stable_id;
    next if $genes{$gene_stable_id} ++; # Process each gene once
    $gene = ( $gene->transform('fpc_pseudomolecule') 
              || $gene->transform('toplevel') || $gene );
    my $gene_seq_region_name = $gene->seq_region_name;
    my $gene_seq_region_start = $gene->seq_region_start;

    $record{gene_stable_id} = $gene_stable_id;
    $record{gene_location} = sprintf
        ("%s:%s", $gene_seq_region_name, $gene_seq_region_start );

    my $gmember = $MEMBER_ADAPTOR->fetch_by_source_stable_id
        ( 'ENSEMBLGENE', $gene_stable_id )
        || die( "[*DIE] Cannot fetch member ENSEMBLGENE $gene_stable_id" );
    my $pmember = $gmember->get_longest_peptide_Member;
    my $pepaligns = get_best_peptide_align_features_by_Member($pmember);

    my( $hpmember );
    my @ogenes;
    my $distance = 9e99;
    my $pepalign;
    foreach my $this_pepalign( @$pepaligns ){
      my $hmember_id = $this_pepalign->{hmember_id};
      my $this_hpmember = &get_hmember_from_pepalign( $this_pepalign );
      my $this_hgmember = $this_hpmember->gene_member;
      my $this_gene = $this_hpmember->get_Gene;
      $this_gene = ( $this_gene->transform('fpc_pseudomolecule') 
              || $this_gene->transform('toplevel') || $gene  );
      my $this_stable_id = $this_gene->stable_id;
      my $this_seq_region_name = $this_gene->seq_region_name;
      my $this_seq_region_start = $this_gene->seq_region_start;
      my $this_location = sprintf
          ("%s:%s", $this_seq_region_name, $this_seq_region_start );

      if( ($this_seq_region_name eq $gene_seq_region_name) ){
        my $this_distance = abs($this_seq_region_start-$gene_seq_region_start);
        if( $this_distance < $distance ){
          # Closest homolog so far
          $record{homolog_stable_id} = $this_stable_id;
          $record{homolog_location} = $this_location;
          $record{distance_to_homolog_bp} = $this_distance;
          $pepalign = $this_pepalign;
          $distance = $this_distance;
          $hpmember = $this_hpmember;
        }
      }
      $record{homolog_stable_id} ||= $this_stable_id;
      $record{homolog_location} ||= $this_location;
      $record{distance_to_homolog_bp} ||= $distance;
      $pepalign ||= $this_pepalign;
      $hpmember ||= $this_hpmember;
      $record{all_homologs} = join
          (', ',($record{all_homologs}||()),"$this_stable_id($this_location)");
    }
    
    $record{gene_length_aa}      = $pmember->seq_length;
    if( $pepalign ){
      my $hgmember = $hpmember->gene_member;
      my $hgene    = $hgmember->get_Gene;
      $record{homolog_length_aa}   = $hpmember->seq_length;
      $record{alignment_length_aa} = $pepalign->{align_length};
      $record{alignment_missmatches} 
      = ( $pepalign->{align_length} - $pepalign->{identical_matches} );
      ($record{alignment_overhang}) 
          = sort{ $b<=>$a } ( ($pepalign->{qstart}-1), ($pepalign->{hstart}-1),
                              ($record{gene_length_aa}-$pepalign->{qend}),
                              ($record{homolog_length_aa}-$pepalign->{hend}) );
    }

    my $tree = $tree_adaptor->fetch_by_Member_root_id($gmember,0);
    unless( $tree ){ # No tree for this gene
      printout( \%record );
      next;
    }
    my $node = $tree->get_leaf_by_Member($gmember)
        || die( "[*DIE] No leaf for ENSEMBLGENE "
                , $gmember->stable_id );

    my $cluster_node = &get_species_specific_cluster_node( $node );
    
    my $gene_count = 0;
    foreach my $leaf( @{$cluster_node->get_all_leaves} ) {
      unless( $leaf->can('gene_member') ){ # Corrupt tree
        next;
      }
      $gene_count++;
    }
    
    $record{tree_node_id} = $tree->node_id;
    $record{tree_leaf_count} = $tree->num_leaves;
    $record{tree_taxon_name} = $tree->get_tagvalue('taxon_name');
    $record{cluster_node_id} = $cluster_node->node_id;
    $record{species_specific_leaf_count} = $gene_count;

    printout( \%record );
    #warn Dumper( \%record );
    $tree->release_tree; # release memory
  }
  #$gene_adaptor->clear_cache; # release memory
  }
}

close OUT;

exit;

#======================================================================
sub printout{
  my $record = shift;
  print OUT join
      ( "\t", 
        map{ defined($record->{$_}) ? $record->{$_} : '' } @FIELDS ) . "\n";
  return 1;
}

#----------------------------------------------------------------------
our $STH;
our $SQL;
sub get_best_peptide_align_features_by_Member{
  my $Member  = shift;
  $SQL ||= qq(
SELECT qmember_id, qstart, qend, hmember_id,hstart, hend, 
       score, evalue, align_length, 
       identical_matches, perc_ident, positive_matches, perc_pos,
       hit_rank, cigar_line
FROM   peptide_align_feature_zea_mays_4577
WHERE  qmember_id = ?
AND    qgenome_db_id = hgenome_db_id
AND    qmember_id != hmember_id );
  $STH ||= $CMP_DBA->dbc->prepare( $SQL );
  my $rv  = $STH->execute( $Member->dbID ) || die $STH->errstr;
  my %ret;
  while( my $rowref = $STH->fetchrow_hashref() ){
    $ret{$rowref->{hit_rank}} ||= [];
    push @{$ret{$rowref->{hit_rank}}}, $rowref;
  }
  my( $best_rank ) = sort{$a<=>$b} keys( %ret );

  return $ret{$best_rank||''} || [];
}
#----------------------------------------------------------------------
sub get_hmember_from_pepalign{
  my $pepalign_hashref = shift || die( "Need a pepalign hashref");
  my $hmember_id = $pepalign_hashref->{hmember_id} ||die("Missing hmember_id");
  my $hmember = $MEMBER_ADAPTOR->fetch_by_dbID($hmember_id);
  return $hmember;
}
#----------------------------------------------------------------------
sub get_peptide_align_feature{
  my $CMP_DBA = shift;
  my $qmember_id = shift;
  my $hmember_id = shift;
  
  my $sql = qq(
SELECT qmember_id, qstart, qend, hmember_id, hstart, hend, 
       score, evalue, align_length, 
       identical_matches, perc_ident, positive_matches, perc_pos,
       hit_rank, cigar_line
FROM   peptide_align_feature_zea_mays_4577
WHERE  qmember_id = ? and hmember_id = ? );

  my $sth = $CMP_DBA->dbc->prepare( $sql );
  my $rv  = $sth->execute( $qmember_id, $hmember_id ) || die $sth->errstr;
  
  return $sth->fetchrow_hashref;
}

#----------------------------------------------------------------------
# Takes a ProteinTree node, and recursively descends the tree until
# other species are seen, or the root node is reached. Returns the
# lowest node that contains all the same species
sub get_species_specific_cluster_node{
  my $node    = shift || return;
  my @leaves = @{$node->get_all_leaves};  
  my $species = $leaves[0]->genome_db->name;
  if( grep{ $_->can('genome_db') and ($_->genome_db->name ne $species) } 
      @leaves ){
    return ; # Mixed species
  }
  return get_species_specific_cluster_node( $node->parent ) || $node;
}


1;


#  LocalWords:  pept
