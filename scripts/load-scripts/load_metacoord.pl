#!/lab/bin/perl -w

BEGIN {
    $ENV{'GrameneDir'} ||= '/usr/local/apache/'; 
    $ENV{'GrameneEnsemblDir'} ||= '/usr/local/apache/'; 
}

# The first shall be last...
use lib map { $ENV{'GrameneDir'}."/$_" } qw ( lib/perl );
use lib map { $ENV{'GrameneEnsemblDir'}."/ensembl-live/$_" } 
        qw ( bioperl-live modules ensembl/modules ensembl-external/modules
             ensembl-draw/modules ensembl-compara/modules );


use strict;
use warnings;


use Getopt::Long;
use Pod::Usage;

use Gramene::Config;
use DBI;

use Data::Dumper qw(Dumper); # For debug

use FindBin qw( $Bin );
use File::Basename qw( dirname );

 
# Import EnsEMBL modules
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Analysis;

=head1 NAME

load_metacoord.pl -- load meta_coord info for the specified species
 and feature table, if a record already exist for the same
 coord_system and feature table, update the max_length, otherwise
 create a new record

=head1 SYNOPSIS

load_metacoord.pl  [options] 
 
 Options:
    --species		species in EnsEMBL registry to use for db [required]
    --v 		makes more verbose
    --help		help message
    --man		full documentation
    --coord_system	coord system to use
    --registry_file	Default is $GrameneEnsemblDir/conf/ensembl.registry
    --table	        name of the ensembl core table that you want to set the coord_system meta info for
    

=head1 OPTIONS

=over 4

=item B<--species> 

A species in the EnsEMBL registry

=item B<--coord_system>

Coordinate system name (e.g. 'clone', 'chromosome' )

Defaults to the sequence level coordinate system 

Can be over-ridden by data file: if region name is chr_NUMBER
perhaps followed by _NUMBER then chromosome coordinates are used
and the 1st NUMBER is the chromosome name.

=item B<--registry_file>

Default is $GrameneEnsemblDir/conf/ensembl.registry

=item B<--table>

name of the ensembl core table that you want to set the coord_system meta info for, can be repeated,
usually

| density_feature
| dna_align_feature     |
| exon                  |
| gene                  |
| karyotype             |
| marker_feature        |
| misc_feature          |
| prediction_exon       |
| prediction_transcript |
| repeat_feature        |
| transcript   



=item B<--help> 

print a help message and exit

=item B<--man> 

print documentation and exit

    
=back



=cut

use vars qw($ENS_DBA $species @tables $coord_system $verbose);
$verbose=0;

{  #Argument Processing
  my $help=0;
  my $man=0;
  my $registry_file;
  GetOptions
      ( 
        "help|?"          => \$help,
        "man"             => \$man,
        "species=s"       => \$species,
	"coord_system=s"  => \$coord_system,
        "registry_file=s" => \$registry_file,

        "table=s"         => \@tables,

        "v+"		  => \$verbose,
        )
      or pod2usage(2);
  pod2usage(-verbose => 2) if $man;
  pod2usage(1) if $help;

  # Validate the input files
  $registry_file    ||= $ENV{GrameneEnsemblDir}.'/conf/ensembl.registry';

  # Load the ensembl file
  $species || ( warn( "Need a --species\n" ) && pod2usage(1) );
  Bio::EnsEMBL::Registry->load_all( $registry_file );
  $ENS_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'core' );
  $ENS_DBA || ( warn( "No core DB for $species set in $registry_file\n" ) &&
                pod2usage(1) );


}

my $csid_sql = q[select coord_system_id from coord_system where name = '%s'];
my $exist_sql = q[select * from meta_coord where coord_system_id=%d and table_name='%s'];
my $max_length_sql = q[select max(seq_region_end-seq_region_start) from %s];
my $meta_coord_update_sql = q[
update  meta_coord 
set     max_length = ? 
where   table_name = ? 
and     coord_system_id = ?
];
my $meta_coord_insert_sql = q[
insert ignore into meta_coord
(table_name, coord_system_id, max_length) 
values (?,?,?) 
];

my $dbh = $ENS_DBA->dbc->db_handle;
$dbh->{RaiseError} = 1;
my $update_sth; 
my $insert_sth; 

my ($csid) = $dbh->selectrow_array( sprintf $csid_sql, $coord_system);
die "Failed select $csid_sql on $coord_system, $DBI::errstr" unless $csid;

for my $tb(sort @tables){

  my ($max_len) = $dbh->selectrow_array( sprintf $max_length_sql, $tb);
  unless( defined $max_len){
    print STDERR "Failed $max_length_sql on $tb, $DBI::errstr";
    next;
  }

  my @result=$dbh->selectrow_array( sprintf $exist_sql, $csid, $tb);
  if (@result){
    unless ($update_sth){
      $update_sth = $dbh->prepare($meta_coord_update_sql) 
	or die "cannot prepare $meta_coord_update_sql, $DBI::errstr" ;
    }

    $update_sth->execute($max_len, $tb, $csid);
    print "succeeded $meta_coord_update_sql on $max_len, $tb, $csid\n";
    
  }else{

    unless ($insert_sth){
      $insert_sth = $dbh->prepare($meta_coord_insert_sql) 
	or die "cannot prepare $meta_coord_insert_sql, $DBI::errstr" ;
    }

    $insert_sth->execute($tb, $csid, $max_len);
    print "succeeded $meta_coord_insert_sql on $tb, $csid, $max_len\n";
  }
  
  
}
