#!/bin/env/perl -w

=pod

=head1 NAME

dump_species_homology_matrix.pl - Generates matrix of cross-species homologs

=head1 SYNOPSIS

  perl dump_species_homology_matrix.pl [options] > resultfile.dat

Options:

  -h|--help             Show brief help and exit.
  -m|--man              Show detailed help
  -e|--ensembl_registry Path to Ensembl registry file.

=head1 OPTIONS

B<-h|--help>
  Print a brief help message and exits.

B<-m|--man>
  Print man page and exit

B<-e|-ensembl_registry>
  Use this Ensembl registry file for database connection info.
  Default is <ENSEMBLHOME>/conf/ensembl.registry

=head1 DESCRIPTION

Script that dumps a pretty-print matrix of cross-species holologs
based on gene tree data. Takes a while to complete...

B<The Ensembl Registry>

  The database connection details for both Ensembl species and compara
  databases must be configured in an ensembl registry file.

  An example registry file would be;
  ---
  use Bio::EnsEMBL::DBSQL::DBAdaptor; 
  use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
  use Bio::EnsEMBL::Registry;
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
  $BASEDIR = dirname(dirname(dirname($Bin)));
  -d $BASEDIR.'/ensembl-live'
    || die( "\n[*DIE] Need $BASEDIR/ensembl-live symlinked to ensembl\n" );
  unshift @INC, $BASEDIR.'/ensembl-live/ensembl/modules';
  unshift @INC, $BASEDIR.'/ensembl-live/ensembl-compara/modules';
}

use Bio::EnsEMBL::Registry;

our $CMP_DBA;
our $CMP_DBH;
BEGIN{  #Argument Processing
  my $help=0;
  my $man=0;
  my( $species, $reg, $slice_name );
  GetOptions
      ( 
        "help|?"             => \$help,
        "man"                => \$man,
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

  $CMP_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( 'compara','compara')
      || pod2usage("\n[*DIE] No compara DB set in $reg\n" );
  $CMP_DBH = $CMP_DBA->dbc;
}

MAIN:{

  # Get a lits of all species plus 
  my %gene_counts = &get_gene_counts_by_species;
  warn Dumper( \%gene_counts );
  my @sp_list = sort keys %gene_counts;

  # Get the homology data for each species
  my %homology_data;
  foreach my $species( @sp_list ){
    my %this_homology_data = &get_homology_summary_for_species($species);
    $homology_data{$species} = {%this_homology_data};
    warn "[INFO] processed homologs for $species\n";
  }
  # Get the orthology data for each species
  my %orthology_data;
  foreach my $species( @sp_list ){
    my %this_orthology_data = &get_orthology_summary_for_species($species);
    $orthology_data{$species} = {%this_orthology_data};
    warn "[INFO] processed orthologs for $species\n";
  }

#  warn Dumper( \%homology_data );

  # Work out the widths of the columns
  foreach my $type ('HOMOLOGY','ORTHOLOGY'){
    my %_data;
    %_data = %homology_data  if $type eq 'HOMOLOGY';
    %_data = %orthology_data if $type eq 'ORTHOLOGY';
    my @colwidths;
    my $ROW = 2;
    my $minwidth = 3;
    foreach my $sp( @sp_list ){
      my( $colwidth ) = ( sort{ $b <=> $a }
                          map{ length($_data{$_}->{$sp}->[$ROW] || 0) }
                          @sp_list );
      $colwidth = 3 if $colwidth < 3;
      push( @colwidths, $colwidth );
    }
    
    # Create a printf template to keep the row spacing constant
    my $rowtmpl = "%3s "
        . join(' ', ( map{"%${_}d"} @colwidths ) ) 
        . "\n";
    
    # Print the header
    my $headtmpl = $rowtmpl;
    $headtmpl =~ tr/d/s/;
    printf( "%s COUNTS\n", uc($type) );
    printf( $headtmpl, '   ', map{&species_abbrev($_)} @sp_list );
    
    # Print the matrix
    foreach my $species( @sp_list ){ 
      my $sp = &species_abbrev($species);
      printf( $rowtmpl, 
              ( $sp, 
                ( map{ $_data{$species}->{$_}->[$ROW] ||0 } 
                  @sp_list ) ) );
    }
  }
}

exit;

#======================================================================
# Turn a species binomial into a 3letter code
sub species_abbrev{
  my $species = shift;
  my $sp;
  
  if( $species =~ /^(.)\S+[\s_](.)\S+[\s_](.)/ ){ # E.g. Oryza sativa Indica -> OsI
    $sp = $1.$2.$3;
  }
  elsif( $species =~ /^(.)\S+[\s_](..)/ ){ # E..g. Zea mays => Zma
    $sp = $1.$2;
  }
  else{ $sp = substr( $species, 0, 3 ) } # Dunno - just take first 3 letters.
  return $sp;
}

#----------------------------------------------------------------------
sub get_gene_counts_by_species{

  my $sql = qq(
select db.name, db.genome_db_id, m.source_name, count(*)
from   genome_db db join member m using( genome_db_id )
where  m.source_name = 'ENSEMBLGENE'
group by 1, 2, 3  );

  my $sth = $CMP_DBH->prepare( $sql );
  my $rv = $sth->execute || die( $sth->errstr );
  my %gene_counts_by_species;
  while( my $row = $sth->fetchrow_arrayref ){
    $gene_counts_by_species{$row->[0]} = $row->[3]; 
  }
  return %gene_counts_by_species;
}

#----------------------------------------------------------------------
sub get_homology_summary_for_species{

  my $species = shift || die( "Need a species name" );

  my $sql = qq(
select this_db.name
     , other_db.name
     , count( distinct( this_m.member_id ) ) number
#     , avg( this_hm.perc_id) this_percid
#     , avg( this_hm.perc_cov) this_perccov
#     , avg( this_hm.perc_id) other_percid
#     , avg( this_hm.perc_cov) other_perccov
from   homology h,
       homology_member this_hm,
       homology_member other_hm,
       member this_m,
       member other_m,
       genome_db this_db,
       genome_db other_db
where h.homology_id=this_hm.homology_id
and   this_hm.member_id=this_m.member_id
and   this_m.genome_db_id=this_db.genome_db_id
and   h.homology_id=other_hm.homology_id
and   other_hm.member_id=other_m.member_id
and   other_m.genome_db_id=other_db.genome_db_id
and   this_hm.member_id != other_hm.member_id
and   this_db.name=?
group by other_db.name
order by other_db.name
);

  my $sth = $CMP_DBH->prepare( $sql );
  my $rv = $sth->execute($species) || die( $sth->errstr );
  my %homology_data;
  while( my @row = $sth->fetchrow_array ){
    $homology_data{$row[1]} = [@row];
  }
  return %homology_data;
}

#----------------------------------------------------------------------
sub get_orthology_summary_for_species{

  my $species = shift || die( "Need a species name" );

  my $sql = qq(
select this_db.name
     , other_db.name
     , count( distinct( this_m.member_id ) ) number
#     , avg( this_hm.perc_id) this_percid
#     , avg( this_hm.perc_cov) this_perccov
#     , avg( this_hm.perc_id) other_percid
#     , avg( this_hm.perc_cov) other_perccov
from   homology h,
       homology_member this_hm,
       homology_member other_hm,
       member this_m,
       member other_m,
       genome_db this_db,
       genome_db other_db
where h.homology_id=this_hm.homology_id
and   this_hm.member_id=this_m.member_id
and   this_m.genome_db_id=this_db.genome_db_id
and   h.homology_id=other_hm.homology_id
and   other_hm.member_id=other_m.member_id
and   other_m.genome_db_id=other_db.genome_db_id
and   this_hm.member_id != other_hm.member_id
and   this_db.name=?
and   h.description in('ortholog_one2one','apparent_ortholog_one2one','ortholog_one2many','ortholog_many2many','possible_ortholog')
group by other_db.name
order by other_db.name
);

  my $sth = $CMP_DBH->prepare( $sql );
  my $rv = $sth->execute($species) || die( $sth->errstr );
  my %homology_data;
  while( my @row = $sth->fetchrow_array ){
    $homology_data{$row[1]} = [@row];
  }
  return %homology_data;
}

#======================================================================
1;
