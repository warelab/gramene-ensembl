#!/usr/bin/env perl

=pod

=head1 NAME

check_genetree_data.pl - QCs the Compara GeneTree data

=head1 SYNOPSIS

  perl check_genetree_data.pl [options]

Options:

 -h|--help  Show brief help and exit.
 -m|--man   Show detailed help
 -u|--url   URL-style connection params to compara DB.
 -l|--long  Run extended test suite.

=head1 OPTIONS

B<-h|--help>
  Print a brief help message and exits.

B<-m|--man>
  Print man page and exit

B<-u|--url>
  URL-style connection params to compara DB in following format:
  mysql://<user>:<pass>@<host>:<port>/<db_name>

B<-l|--long>
  Run extended test suite.

=head1 DESCRIPTION

  Add a description of each test here.

Maintained by Albert Vilella <avilella@ebi.ac.uk>

=cut

use strict;
use warnings;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Hive::URLFactory;

use Getopt::Long;
use Pod::Usage;

my $DEFAULT_URL = 'mysql://ensro@127.0.0.1:4313/mm14_compara_homology_68';

# Get options
my $help=0;
my $man=0;
my( $url, $long, $V );
  GetOptions
      ( 
        "help|?"             => \$help,
        "man"                => \$man,
        "url=s"              => \$url,
        "longtests=s"        => \$long,
        "verbose"            => \$V, # Not yet used
        )
    or pod2usage(2);
pod2usage(-verbose => 2) if $man;
pod2usage(1) if $help;

$url ||= $DEFAULT_URL;

my $dba = Bio::EnsEMBL::Hive::URLFactory->fetch($url . ';type=compara');

my $doit = 0;
$doit = 1 if ($long);

$|=1;

my ($sql, $sth);

if ($doit) {

# # Check data consistency between gene_count and number of homology entries
# ##########################################################################

#   print "# Check data consistency between gene_count and number of homology entries\n";

#   $sql = "select node_id,value,(value*(value-1))/2 from gene_tree_tag where tag='gene_count'";

#   $sth = $dba->dbc->prepare($sql);
#   $sth->execute;

#   my $this_count = sprintf("%05d",1);
#   while (my $aref = $sth->fetchrow_arrayref) {
#     my ($node_id, $value, $homology_count) = @$aref;
#     my $sql2 = "select count(*) from homology where tree_node_id=$node_id";
#     my $sth2 = $dba->dbc->prepare($sql2);
#     print STDERR "## Checking $node_id / gene_count = $value / $homology_count / [$this_count]\n";
#     $this_count = sprintf("%05d",$this_count+1);
#     $sth2->execute;
#     my $count = $sth2->fetchrow_array;
#     if ($count !=0 && $count != $homology_count) {
#       print STDERR "ERROR: tree $node_id (gene_count = $value) gene_count != homologies : should have $homology_count homologies instead of $count\n";
#       print STDERR "ERROR: USED SQL : $sql\n                  $sql2\n";
#     }
#     $sth2->finish;

#     my $sql3 = "select count(*) from homology h, homology_member hm where h.homology_id=hm.homology_id and h.tree_node_id=$node_id";
#     my $sth3 = $dba->dbc->prepare($sql3);
#     $sth3->execute;
#     $count = $sth3->fetchrow_array;
#     if ($count !=0 && $count != 2*$homology_count) {
#       print STDERR "ERROR: tree $node_id (gene_count = $value) gene_count != homology_member : should have $homology_count homologies instead of $count\n";
#       print STDERR "ERROR: USED SQL : $sql\n                  $sql3\n";
#     }
#     $sth3->finish;
#   }

#   $sth->finish;

# Check for dangling internal nodes that have no children
######################################################
print "# Check for dangling internal nodes that have no children\n";

  $sql = "select count(*) from gene_tree_node n1 left join gene_tree_node n2 on n1.node_id=n2.parent_id where n2.parent_id is NULL and n1.right_index-n1.left_index > 1";
  
  $sth = $dba->dbc->prepare($sql);
  $sth->execute;
  
  while (my $aref = $sth->fetchrow_arrayref) {
    #should 0, if not delete culprit node_id in gene_tree_member
    my ($count) = @$aref;
    if ($count == 0) {
      print "PASSED: gene_tree_node is consistent - no dangling internal nodes\n";
    } else {
      print STDERR "ERROR: gene_tree_node has dangling internal nodes with no children based on the left and right_index\n";
      print STDERR "ERROR: USED SQL : $sql\n";
    }
  }
  
  $sth->finish;

## This check is not relevant any more since we don't store all between_species_paralogs
# # Check data consistency between pt* tables on node_id
# ######################################################
# print "# Check data consistency between pt* tables on node_id\n";

#   $sql = "select count(*) from gene_tree_member ptm left join gene_tree_node ptn on ptm.node_id=ptn.node_id where ptn.node_id is NULL";
  
#   $sth = $dba->dbc->prepare($sql);
#   $sth->execute;
  
#   while (my $aref = $sth->fetchrow_arrayref) {
#     #should 0, if not delete culprit node_id in gene_tree_member
#     my ($count) = @$aref;
#     if ($count == 0) {
#       print "PASSED: gene_tree_member versus gene_tree_node is consistent\n";
#     } else {
#       print STDERR "ERROR: gene_tree_member versus gene_tree_node is NOT consistent\n";
#       print STDERR "ERROR: USED SQL : $sql\n";
#     }
#   }
  
#   $sth->finish;
  
#   $sql = "select count(*) from gene_tree_tag ptt left join gene_tree_node ptn on ptt.node_id=ptn.node_id where ptn.node_id is NULL";
  
#   $sth = $dba->dbc->prepare($sql);
#   $sth->execute;
  
#   while (my $aref = $sth->fetchrow_arrayref) {
#     #should be 0, if not delete culprit node_id in gene_tree_tag
#     my ($count) = @$aref;
#     if ($count == 0) {
#       print "PASSED: gene_tree_tag versus gene_tree_node is consistent\n";
#     } else {
#       print STDERR "ERROR: gene_tree_tag versus gene_tree_node is NOT consistent\n";
#       print STDERR "ERROR: USED SQL : $sql\n";
#     }
#   }
  
#   $sth->finish;

# # Check that the gene_count tags and the right-left index counts are equivalent
#   $sql = "select p1.*, ROUND((((p1.right_index-p1.left_index+1)/2)+1)/2) as p1_size, ptt1.value as gene_count from gene_tree_node p1, gene_tree_tag ptt1 where p1.parent_id=1 and p1.node_id=ptt1.node_id and ptt1.tag='gene_count' and ROUND((((p1.right_index-p1.left_index+1)/2)+1)/2)!=ptt1.value";
  
#   $sth = $dba->dbc->prepare($sql);
#   $sth->execute;
  
#   while (my $aref = $sth->fetchrow_arrayref) {
#     #should be 0, if not delete culprit node_id in gene_tree_tag
#     my ($count) = @$aref;
#     if ($count == 0) {
#       print "PASSED: gene_tree_tag gene_count versus right_index left_index formula is consistent\n";
#     } else {
#       print STDERR "ERROR: gene_tree_tag versus gene_tree_node is NOT consistent\n";
#       print STDERR "ERROR: USED SQL : $sql\n";
#     }
#   }
  
#   $sth->finish;


# Check for unique member presence in ptm
#########################################
print "# Check for unique member presence in ptm\n";

  #$sql = "select member_id from gene_tree_member group by member_id having count(*)>1";
  $sql = "select member_id from homology_member group by member_id having count(*)>1";

  #This should return 0 rows. 
  
  $sth = $dba->dbc->prepare($sql);
  
  # If no duplicate, each select should return an empty row;
  $sth->execute;
  
  my $ok = 1;
  while (my $aref = $sth->fetchrow_arrayref) {
    my ($member_id) = @$aref;
    print STDERR "ERROR: some member duplicates in gene_tree_member! This can happen if there are remains of broken clusters in the db that you haven't deleted yet. Usually before merging to the compara production db\n";
    print STDERR "ERROR: USED SQL : $sql\n";
    $ok = 0;
    last;
  }
  print "PASSED: no member duplicates in gene_tree_member\n" if ($ok);
  

## # check for unique member presence in ptm but without considering singletons
## ############################################################################
## 
##   $sql = 'select * from gene_tree_member ptm, gene_tree_node ptn where 0<abs(ptn.left_index-ptn.right_index) and ptn.node_id=ptm.node_id group by ptm.member_id having count(*)>1';
##   
##   #This should return 0 rows. 
##   
##   $sth = $dba->dbc->prepare($sql);
##   
##   # If no duplicate, each select should return an empty row;
##   $sth->execute;
##   
##   $ok = 1;
##   while (my $aref = $sth->fetchrow_arrayref) {
##     my ($member_id) = @$aref;
##     print STDERR "ERROR: some member duplicates (non-singletons) in gene_tree_member!\n";
##     print STDERR "ERROR: USED SQL : $sql\n";
##     $ok = 0;
##     last;
##   }
##   print "PASSED: no member duplicates (non-singletons) in gene_tree_member\n" if ($ok);
##   
## 

# Check data consistency between pt_node and homology with node_id
##################################################################
print "# Check data consistency between pt_node and homology with node_id\n";


  #$sql = "select count(*) from homology h left join gene_tree_node ptn on h.ancestor_node_id=ptn.node_id where h.description!='other_paralog' AND ptn.node_id is NULL";
  
   $sql = "select count(*) from homology h left join gene_tree_node ptn on h.gene_tree_node_id=ptn.node_id where h.description!='other_paralog' AND ptn.node_id is NULL";

  $sth = $dba->dbc->prepare($sql);
  $sth->execute;
  
  while (my $aref = $sth->fetchrow_arrayref) {
    #should be 0
    my ($count) = @$aref;
    if ($count == 0) {
      print "PASSED: homology versus gene_tree_node is consistent\n";
    } else {
      print STDERR "ERROR: homology versus gene_tree_node is NOT consistent\n";
      print STDERR "ERROR: USED SQL : $sql\n";
    }
  }
  
  $sth->finish;

} # end of if ($doit)

# Check for homology has no duplicates
######################################
print "# Check for homology has no duplicates\n";

my $mlsses;
my $mlssa;

$mlssa = $dba->get_MethodLinkSpeciesSetAdaptor;
$mlsses = $mlssa->fetch_all_by_method_link_type('ENSEMBL_ORTHOLOGUES');

push @{$mlsses},@{$mlssa->fetch_all_by_method_link_type('ENSEMBL_PARALOGUES')};

$sql = "select hm1.member_id,hm2.member_id,h.method_link_species_set_id from homology_member hm1, homology_member hm2, homology h where h.homology_id=hm1.homology_id and hm1.homology_id=hm2.homology_id and hm1.member_id<hm2.member_id and h.method_link_species_set_id=? group by hm1.member_id,hm2.member_id having count(*)>1";

$sth = $dba->dbc->prepare($sql);

while (my $mlss = shift @$mlsses) {
  # If no duplicate, each select should return an empty row;
  $sth->execute($mlss->dbID);
  my $ok = 1;
  while (my $aref = $sth->fetchrow_arrayref) {
    my ($member_id1, $member_id2,$mlss_id) = @$aref;
    print STDERR "ERROR: some homology duplicates in method_link_species_set_id=$mlss_id!\n";
    print STDERR "ERROR: USED SQL : $sql\n";
    $ok = 0;
    last;
  }
  print "PASSED: no homology duplicates in method_link_species_set_id=".$mlss->dbID."\n" if ($ok);
}

$sth->finish;

# Check that one2one genes are not involved in any other orthology
# (one2many or many2many)
# check to be done by method_link_species_set with method_link_type
# ENSEMBL_ORTHOLOGUES
##################################################################
print "# Check that one2one genes are not involved in any other orthology\n";

$mlssa = $dba->get_MethodLinkSpeciesSetAdaptor;
$mlsses = $mlssa->fetch_all_by_method_link_type('ENSEMBL_ORTHOLOGUES');

if ($doit) {

# $sql = "select h.method_link_species_set_id,hm.member_id from homology h, homology_member hm where h.method_link_species_set_id=? and h.homology_id=hm.homology_id group by hm.member_id having count(*)>1 and group_concat(h.description) like '%ortholog_one2one%'";
$sql = "select h.method_link_species_set_id,hm.member_id from homology h, homology_member hm where h.method_link_species_set_id=? and h.homology_id=hm.homology_id group by hm.member_id having count(*)>1 and group_concat(h.description) ='ortholog_one2one'";

$sth = $dba->dbc->prepare($sql);

while (my $mlss = shift @$mlsses) {
  # If no mistake in OrthoTree, each select should return an empty row;
  $sth->execute($mlss->dbID);
  my $ok = 1;
  while (my $aref = $sth->fetchrow_arrayref) {
    my ($mlss_id, $member_id) = @$aref;
    # print STDERR "ERROR: some [apparent_]one2one_ortholog also one2many or many2many in method_link_species_set_id=$mlss_id!\n";
    print STDERR "ERROR: some one2one_ortholog also one2many or many2many in method_link_species_set_id=$mlss_id!\n";
    print STDERR "ERROR: USED SQL : $sql\n";
    $ok = 0;
    last;
  }
  # print "PASSED: [apparent_]one2one_ortholog are ok in method_link_species_set_id=".$mlss->dbID."\n" if ($ok);
  print "PASSED: one2one_ortholog are ok in method_link_species_set_id=".$mlss->dbID."\n" if ($ok);
}

$sth->finish;

}
