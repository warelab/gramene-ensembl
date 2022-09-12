#!/usr/bin/perl -w

#retrieve meta data from core database and update the metadatabase

use lib "/usr/local/ensembl-87/ensemblgenomes-api/modules/";
use lib "/usr/local/ensembl-87/ensembl-production/modules/";
use lib "/usr/local/ensembl-live/ensembl-hive/modules/";


use strict;
use Getopt::Long;
use Data::Dumper;
use Pod::Usage;
use File::Basename;
use Cwd 'abs_path';
use vars qw[ $VERSION ];
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::DBSQL::DBConnection; 
use Bio::EnsEMBL::Utils::MetaData::DBSQL::GenomeInfoAdaptor;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Production::Pipeline::Production::GenomeStats;
use Bio::EnsEMBL::Production::Pipeline::Production::StatsGenerator;
use Bio::EnsEMBL::Hive::Process;
use Bio::EnsEMBL::Hive::AnalysisJob;

pod2usage(0) unless $ARGV[0];

my ($_help, $dbhost, $dbport, $dbuser, $dbpass, $cdbname, $mdbname, $write);
my ($reg, $species);

GetOptions(
           'h|help'        => \$_help,
           'dbhost=s'      => \$dbhost,
           'dbport=s'      => \$dbport,
           'dbuser=s'      => \$dbuser,
           'dbpass=s'      => \$dbpass,
           'cdbname=s'     => \$cdbname,
           'mdbname=s'     => \$mdbname,
           'write'          => \$write,
	   'reg=s'        => \$reg,
	   'species=s'    => \$species,    
       ) or die;

pod2usage(2) if $_help;

#pod2usage(2) unless ( $dbhost && $dbuser && $dbpass && $cdbname && $mdbname );
#$dbport ||= '3306';
#print STDERR "Target: $dbhost $dbport $cdbname $mdbname\n";

pod2usage(2) unless ( $reg && $species);

Bio::EnsEMBL::Registry->load_all( $reg );
my  $cdba = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'core' );
  $cdba || ( print( "No core DB for $species set in $reg\n" ) &&
                pod2usage(2) );



## Get core db DBAdaptor
#my $cdba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
#                                               -user   => $dbuser,
#                                               -pass   => $dbpass,
#                                               -dbname => $cdbname,
#                                               -host   => $dbhost,
#                                               -port   => $dbport,
#                                               -driver => 'mysql',
#                                               );

my %genome_info;

$genome_info{name}            = $genome_info{species} = $cdba->get_MetaContainer()->get_production_name;
$genome_info{taxonomy_id}     = $cdba->get_MetaContainer()->get_taxonomy_id();
$genome_info{division}        = $cdba->get_MetaContainer()->get_division() || 'Ensembl';
$genome_info{assembly_name}   = $cdba->get_MetaContainer()->single_value_by_key('assembly.name');
$genome_info{assembly_level}  = 'chromosome';
$genome_info{genebuild}       = $cdba->get_MetaContainer()->single_value_by_key('genebuild.version');
$genome_info{dbname}          = $cdba->dbc()->dbname() ;

#Bio::EnsEMBL::Registry->add_adaptor($genome_info{name}, "core", "GenomeContainer", $cdba);
#my $genome = Bio::EnsEMBL::Registry->get_adaptor( $genome_info{name}, "core", "GenomeContainer" );

my $genome_obj =
    Bio::EnsEMBL::Registry->get_adaptor( $species, "core", "GenomeContainer" );

$genome_info{base_count} = $genome_obj->get_total_length();

unless ( $genome_info{base_count} ) {

 my $stats_obj = new Bio::EnsEMBL::Production::Pipeline::Production::GenomeStats;
    #Bio::EnsEMBL::Production::Pipeline::Production::StatsGenerator;

  &run($stats_obj, $species);
 print "Updated stats for $species\n";

 my $genome_obj2 =
    Bio::EnsEMBL::Registry->get_adaptor( $species, "core", "GenomeContainer" );

    $genome_info{base_count} = $genome_obj2->get_total_length();

   die "Failed to update stats and get base_count" unless ( $genome_info{base_count} );
}

#get_ref_length

for my $k (keys %genome_info) {

	print "$k => $genome_info{$k}\n";

}

my $genome = Bio::EnsEMBL::Utils::MetaData::GenomeInfo->new( 
    -name             => $genome_info{name},
    -species          => $genome_info{species},
    -species_id       => 1,
    -taxonomy_id      => $genome_info{taxonomy_id},
    -species_taxonomy_id => $genome_info{taxonomy_id},
    -division         => $genome_info{division},
    -assembly_name    => $genome_info{assembly_name},
    -assembly_level   => $genome_info{assembly_level},
    -genebuild        => $genome_info{genebuild},
    -dbname           => $genome_info{dbname} );

$genome->base_count( $genome_info{base_count} );

my $metadb_dba = Bio::EnsEMBL::Registry->get_DBAdaptor('meta', 'metadata');

#my $genome = $metadb_dba->fetch_by_species('arabidopsis_thaliana');
#print "test => ". $genome->name()."\n";

$metadb_dba->store($genome) if $write;



sub run {
  my $stats_obj = shift;
  my $species    = shift;
  warn("debug weix: species => $species");
  my $dba        = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'core');
  my $ta         = Bio::EnsEMBL::Registry->get_adaptor($species, 'core', 'transcript');
  my $aa         = Bio::EnsEMBL::Registry->get_adaptor($species, 'core', 'attribute');

  my $has_readthrough = 0;
  my @readthroughs = @{ $aa->fetch_all_by_Transcript(undef, 'readthrough_tra') };
  if (@readthroughs) {
    $has_readthrough = 1;
  }

  my $sum = 0;
  my $count;
  my $alt_count;
  my %slices_hash;
  my %stats_hash;
  my %stats_attrib;
  
  my $total = scalar(@{ Bio::EnsEMBL::Registry->get_adaptor($species, 'core', 'gene')->fetch_all });

   #my %genome_counts = (
   # PredictionTranscript => 'core',
   # StructuralVariation => 'variation',
   #);

  my %attrib_codes = $stats_obj->get_attrib_codes( $has_readthrough );
  my %alt_attrib_codes = &get_alt_attrib_codes($stats_obj, $has_readthrough);
  my $slice_adaptor = $dba->get_SliceAdaptor();

  my $all_slices = $slice_adaptor->fetch_all('toplevel');
  my @all_sorted_slices =
   sort( { $a->coord_system()->rank() <=> $b->coord_system()->rank()
           || $b->seq_region_length() <=> $a->seq_region_length() } @$all_slices) ;
  while (my $slice = shift @all_sorted_slices) {
	$stats_hash{'total_length'} += $slice->length;
    if ($slice->is_reference) {
      $stats_hash{'transcript'} += $ta->count_all_by_Slice($slice);
      $stats_attrib{'transcript'} = 'transcript_cnt';
warn("weix debug before attrib_codes loop");
      foreach my $ref_code (keys %attrib_codes) {
	warn("debug ref_code=$ref_code");
        $count = &get_feature_count($stats_obj, $slice, $ref_code, $attrib_codes{$ref_code}, $species);
        if ($count > 0) {
          &store_attrib($stats_obj, $slice, $count, $ref_code, $species);
        }
        $stats_hash{$ref_code} += $count;
      }
warn("weix debug after attrib_codes loop");


      $count = &get_attrib($stats_obj, $slice, 'noncoding_cnt_%', $species);
      if ($count > 0) {
        store_attrib( $stats_obj, $slice, $count, 'noncoding_cnt', $species);
      }
      $stats_hash{'noncoding_cnt'} += $count;
      $stats_attrib{'noncoding_cnt'} = 'noncoding_cnt';
    } elsif ($slice->seq_region_name =~ /LRG/) {
      next;
    } else {
      my $sa = Bio::EnsEMBL::Registry->get_adaptor($species, 'core', 'slice');
      my $alt_slices = $sa->fetch_by_region_unique($slice->coord_system->name(), $slice->seq_region_name());
      foreach my $alt_slice (@$alt_slices) {
        $stats_hash{'alt_transcript'} += $ta->count_all_by_Slice($alt_slice);
        $stats_attrib{'alt_transcript'} = 'transcript_acnt';
        foreach my $alt_code (keys %alt_attrib_codes) {
          $alt_count = &get_feature_count($stats_obj, $alt_slice, $alt_code, $alt_attrib_codes{$alt_code}, $species);
          if ($alt_count > 0) {
            &store_attrib($stats_obj, $slice, $alt_count, $alt_code, $species);
          }
          $stats_hash{$alt_code} += $alt_count;
        }
        $alt_count = &get_attrib($stats_obj, $slice, 'noncoding_acnt_%', $species);
        if ($alt_count > 0) {
          &store_attrib($stats_obj, $slice, $alt_count, 'noncoding_acnt', $species);
        }
        $stats_hash{'noncoding_acnt'} += $alt_count;
        $stats_attrib{'noncoding_acnt'} = 'noncoding_acnt';
      }
    }
    if ($sum >= $total) {
      last;
    }
  }

  $stats_hash{'ref_length'} = $stats_hash{'total_length'};
  &store_statistics($stats_obj, $species, \%stats_hash, \%stats_attrib);
  #Disconnecting from the registry
  $dba->dbc->disconnect_if_idle();
  $ta->dbc->disconnect_if_idle();
  $aa->dbc->disconnect_if_idle();
}


sub get_alt_attrib_codes {
  my ($self, $has_readthrough) = @_;
  my @alt_attrib_codes = ('coding_acnt', 'pseudogene_acnt', 'noncoding_acnt_s', 'noncoding_acnt_l', 'noncoding_acnt_m');
  if ($has_readthrough) {
    push @alt_attrib_codes, ('coding_racnt', 'pseudogene_racnt', 'noncoding_racnt_s', 'noncoding_racnt_l', 'noncoding_racnt_m');
  }
  my %biotypes;
  foreach my $alt_code (@alt_attrib_codes) {
    my ($group, $subgroup) = $alt_code =~ /(\w+)\_r?acnt_?([a-z]?)/;
    if ($subgroup) { $group = $subgroup . $group; }
    my $biotypes = $self->get_biotype_group($group);
    $biotypes{$alt_code} = $biotypes;
  }
  return %biotypes;
}

sub get_feature_count {
  my ($self, $table, $dbtype, $condition, $species) = @_;
  my $dba = Bio::EnsEMBL::Registry->get_DBAdaptor($species, $dbtype);
  my $count = 0;
  if ($table =~ /Length/) {
    my $slice_adaptor = $dba->get_adaptor('slice');
    my $slices = $slice_adaptor->fetch_all('seqlevel');
    foreach my $slice (@$slices) {
      $count += $slice->length();
    }
  } elsif (defined $dba) {
    my $adaptor = $dba->get_adaptor($table);
    $count = $adaptor->generic_count($condition);
  }
  return $count;
}

sub store_attrib {
  my ($self, $slice, $count, $code, $species) = @_;
  my $aa          = Bio::EnsEMBL::Registry->get_adaptor($species, 'core', 'Attribute');
  my $prod_dba    = Bio::EnsEMBL::Registry->get_DBAdaptor('multi', 'production'); #$self->get_production_DBAdaptor();
  my $prod_helper = $prod_dba->dbc()->sql_helper();
  my $sql = q{
    SELECT name, description
    FROM attrib_type
    WHERE code = ? };
  my ($name, $description) = @{$prod_helper->execute(-SQL => $sql, -PARAMS => [$code])->[0]};
  my $attrib = Bio::EnsEMBL::Attribute->new(
    -NAME        => $name,
    -CODE        => $code,
    -VALUE       => $count,
    -DESCRIPTION => $description
  );
  my @attribs = ($attrib);
  $aa->remove_from_Slice($slice, $code);
  $aa->store_on_Slice($slice, \@attribs);
  $prod_dba->dbc()->disconnect_if_idle();
}

sub get_attrib {
  my ($self, $slice, $code, $species) = @_;
  my $aa          = Bio::EnsEMBL::Registry->get_adaptor($species, 'core', 'Attribute');
  my $attributes = $aa->fetch_all_by_Slice($slice, $code);
  my $count = 0;
  foreach my $attribute (@$attributes) {
    $count += $attribute->value();
  }
  return $count;
}

sub store_statistics {
  my ($self, $species, $stats_hash, $stats_attrib) = @_;
  my $stats;
  my $genome_container = Bio::EnsEMBL::Registry->get_adaptor($species, 'core', 'GenomeContainer');
  my $attrib;
  foreach my $stats (keys %$stats_hash) {
    if ($stats_attrib->{$stats}) {
      $attrib = $stats_attrib->{$stats};
    } else {
      $attrib = $stats;
    }
    $genome_container->store($stats, $stats_hash->{$stats}, $attrib);
  }
}

__END__


## EG
  my $metadata_db = $stats_obj->full_tree->{MULTI}->{databases}->{DATABASE_METADATA};
  my $genome_info_adaptor;
  if ($metadata_db) {
    my $dbc = Bio::EnsEMBL::DBSQL::DBConnection->new(
      -USER   => $dbuser,
      -PASS   => $dbpass,
      -PORT   => $dbport,
      -HOST   => $dbhost,
      -DBNAME => $mdbname
    );
    $genome_info_adaptor = Bio::EnsEMBL::Utils::MetaData::DBSQL::GenomeInfoAdaptor->new(-DBC => $dbc);
  }
##


my $genome =
  Bio::EnsEMBL::Utils::MetaData::GenomeInfo->new(
  -name    => $dba->species(),
  -species    => $dba->species(),
  -species_id => $dba->species_id(),
  -taxonomy_id => $dba->get_MetaContainer()->get_taxonomy_id(),
  -division   => $dba->get_MetaContainer()->get_division() || 'Ensembl',
  -assembly_name   => $dba->get_MetaContainer()->single_value_by_key('assembly.name'),
  -assembly_level   => 'spirit',
  -genebuild => $dba->get_MetaContainer()->single_value_by_key('genebuild.version'),
  -dbname     => $dba->dbc()->dbname() );
  $genome->base_count(1664);
ok( defined $genome, "GenomeInfo exists" );
ok(! defined $genome->dbID(), "GenomeInfo has no ID" );

diag("Storing a new GenomeInfo object");
$gdba->store($genome);

