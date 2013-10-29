#!/bin/env perl -w

# A quick-and-dirty script to identify duplicated genes in the maize
# database based on the paralog assignments and alignments between
# them in compara

# This script has to be run from brie to get the paths correct

use strict;

# Some constants
our $REGFILE = '/usr/local/gramene/conf/ensembl.registry';
our $SPECIES = 'Zea_mays';
our $COMPARA = 'compara';
our $METHOD_LINK = 'ENSEMBL_PARALOGUES';
our $MAX_DISTANCE= 1000000;
our $OUTFILE = 'maize_duplicates.tsv';
our @FIELDS  = qw( gene1 
                   gene2 
                   gene1_length_aa 
                   gene2_length_aa
                   sequence_name 
                   distance_bp
                   align_length
                   missmatches
                   percent_id
                   gene1_start
                   gene1_end
                   gene2_start
                   gene2_end );

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
my ($tax_id) = @{$meta_container->list_value_by_key('species.taxonomy_id')};
$tax_id || die( "\n[*DIE] No meta species.taxonomy_id in core DB!\n" );
my $genomedb_adaptor = $CMP_DBA->get_adaptor('GenomeDB');
my $mlss_adaptor     = $CMP_DBA->get_adaptor('MethodLinkSpeciesSet');
my $homology_adaptor = $CMP_DBA->get_adaptor('Homology');
my $member_adaptor   = $CMP_DBA->get_adaptor('Member');

my( $genome_db ) = ( grep{ $_->taxon_id == $tax_id } 
                     @{$genomedb_adaptor->fetch_all} );
$genome_db || die( "\n[*DIE] No taxon_id $tax_id in compara.\n" );

my $mlss = $mlss_adaptor->fetch_by_method_link_type_GenomeDBs
    ( $METHOD_LINK, [$genome_db] )
    || die( "\n[*DIE] No mlss for $METHOD_LINK and taxon $tax_id.\n" );

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
  
  my @homologs = @{$homology_adaptor
                       ->fetch_all_by_MethodLinkSpeciesSet($mlss)};
  #my @homologs = ($homology_adaptor->fetch_by_dbID(2));

  while( my $homology = shift @homologs  ){
    my %record;

    my @member_pair = @{$homology->gene_list};

    # Retrieve the Ensembl genes fom the homologs, and project to toplevel
    my @clone_gene_pair = ( map{$_->get_Gene} @member_pair );
    my @gene_pair;
    my @pept_pair;
    my $untransformable = 0;
    foreach my $gene( @clone_gene_pair ){
      my $fpc_gene = $gene->transform('fpc_pseudomolecule');
      unless( $fpc_gene ){
#        warn( "[WARN] ", $gene->seq_region_name, 
#              " does not project to fpc_pseudomolecule coord system\n" );
        $untransformable{$gene->stable_id} = $gene->seq_region_name;
        $untransformable++
      }
      push @gene_pair, $fpc_gene;
    }
    next if $untransformable;

    foreach my $member_attribute (@{$homology->get_all_Member_Attribute}) {
      my ($member, $attribute) = @{$member_attribute};
      my $peptide_member = $member_adaptor->fetch_by_dbID
          ($attribute->peptide_member_id);
      push @pept_pair, $peptide_member;
    }

    # Are both genes on the same FPContig? Within threshold BP?
    next if  $gene_pair[0]->seq_region_name ne $gene_pair[1]->seq_region_name;
    $record{distance_bp} = abs( $gene_pair[0]->seq_region_start 
                             - $gene_pair[1]->seq_region_start );
    next if $record{distance_bp} > $MAX_DISTANCE;
    
    $record{gene1}           = $gene_pair[0]->stable_id;
    $record{gene2}           = $gene_pair[1]->stable_id;
    $record{gene1_length_aa} = $pept_pair[0]->seq_length;
    $record{gene2_length_aa} = $pept_pair[1]->seq_length;
    $record{sequence_name}   = $gene_pair[0]->seq_region_name;
    
#    warn( "  > Genes $record{gene1} and $record{gene2} "
#          ," on seq $record{sequence} "
#          , " separated by $record{distance} bp \n");

    # OK - genes are in a location consistent with duplication. 
    # Is the alignment consistent?

    my $pep_align_feat = &get_peptide_align_feature
        ( $CMP_DBA, map{$_->member_id} @pept_pair );
    unless( $pep_align_feat ){
      # Not in top-10-ranked hit. Ignore
#      warn( "[WARN] No peptide_align_feature for member_ids ",
#           join( ', ', map{$_->member_id} @pept_pair ) );
      next;
    }
    
    $record{align_length} = $pep_align_feat->{align_length};
    $record{percent_id}   = $pep_align_feat->{perc_ident};
    $record{gene1_start}  = $pep_align_feat->{qstart};
    $record{gene1_end}    = $pep_align_feat->{qend};
    $record{gene2_start}  = $pep_align_feat->{hstart};
    $record{gene2_end}    = $pep_align_feat->{hend};
    $record{missmatches}  = ( $pep_align_feat->{align_length} 
                              - $pep_align_feat->{identical_matches} );
    $record{gene1_overhang5} = $record{gene1_start} - 1;
    $record{gene2_overhang5} = $record{gene2_start} - 1;
    $record{gene1_overhang3} = $record{gene1_length_aa} - $record{gene1_end};
    $record{gene2_overhang3} = $record{gene2_length_aa} - $record{gene2_end};
    ($record{gene1_overhang}) = sort{$b<=>$a} ( $record{gene1_overhang5},
                                                $record{gene1_overhang3} );
    ($record{gene2_overhang}) = sort{$b<=>$a} ( $record{gene2_overhang5},
                                                $record{gene2_overhang3} );

    print OUT join
        ( "\t", 
          map{ $record{$_} } @FIELDS ) . "\n";
    
#    use Data::Dumper qw (Dumper);
#    warn Dumper( \%record );

    $dup_count++;
    $seen_genes{$record{gene1}} ++;
    $seen_genes{$record{gene2}} ++;

    if( ! $record{missmatches} ){
      if( ! $record{gene1_overhang} and ! $record{gene2_overhang} ){
        $identical_genes{$record{gene1}} ++;
        $identical_genes{$record{gene2}} ++;
      }
      else{
        $overhanging_genes{$record{gene1}} ++;
        $overhanging_genes{$record{gene2}} ++;
      }
    }
    elsif( $record{missmatches} ==1 ){
      $single_missmatch_genes{$record{gene1}} ++;
      $single_missmatch_genes{$record{gene2}} ++;
    }

    warn( "[INFO] Dup count: $dup_count" );
    
  }
} 

warn "[INFO] Total gene duplicates: " . scalar(keys %seen_genes)."\n";
my %counts = ();
foreach my $gene( keys %seen_genes ){
  $counts{$seen_genes{$gene}} ++;
}
foreach my $count( sort{$a<=>$b} keys %counts ){
  warn "[INFO] ...with $count pairs: $counts{$count}\n";
}

warn "[INFO] Total identical genes: " . scalar(keys %identical_genes)."\n";
%counts = ();
foreach my $gene( keys %identical_genes ){
  $counts{$identical_genes{$gene}} ++;
}
foreach my $count( sort{$a<=>$b} keys %counts ){
  warn "[INFO] ...with $count pairs: $counts{$count}\n";
}

warn "[INFO] Total overhanging genes: " . scalar(keys %overhanging_genes)."\n";
%counts = ();
foreach my $gene( keys %overhanging_genes ){
  $counts{$overhanging_genes{$gene}} ++;
}
foreach my $count( sort{$a<=>$b} keys %counts ){
  warn "[INFO] ...with $count pairs: $counts{$count}\n";
}

warn "[INFO] Total single_missmatch genes: " 
     . scalar(keys %single_missmatch_genes)."\n";
%counts = ();
foreach my $gene( keys %single_missmatch_genes ){
  $counts{$single_missmatch_genes{$gene}} ++;
}
foreach my $count( sort{$a<=>$b} keys %counts ){
  warn "[INFO] ...with $count pairs: $counts{$count}\n";
}

my %unt_clones=map{$untransformable{$_} = $_} keys %untransformable;
warn "[INFO] Could not project to chromosome: " 
    . scalar(keys %untransformable) . " genes from "
    . scalar(keys %unt_clones) . " clones\n";

close OUT;

exit;

#======================================================================
sub get_peptide_align_feature{
  my $CMP_DBA = shift;
  my $qmember_id = shift;
  my $hmember_id = shift;
  
  my $sql = qq(
SELECT qstart, qend, hstart, hend, 
       score, evalue, align_length, 
       identical_matches, perc_ident, positive_matches, perc_pos,
       hit_rank, cigar_line
FROM   peptide_align_feature_zea_mays_4577
WHERE  qmember_id = ? and hmember_id = ? );

  my $sth = $CMP_DBA->dbc->prepare( $sql );
  my $rv  = $sth->execute( $qmember_id, $hmember_id ) || die $sth->errstr;
  
  return $sth->fetchrow_hashref;
}


1;


#  LocalWords:  pept
