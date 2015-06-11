package Reg;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Production::DBSQL::DBAdaptor;

Bio::EnsEMBL::Registry->no_version_check(1);
Bio::EnsEMBL::Registry->no_cache_warnings(1);
{
  my $version = 79;
  Bio::EnsEMBL::Registry->load_registry_from_multiple_dbs(
    {
      -host => "cabot",
      -port => 3306,
      -db_version => $version,
      -user => "weix",
      -pass => 'warelab',
      -NO_CACHE => 1,
    },
  );
  Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -HOST => 'cabot',
    -PORT => 3306,
    -USER => 'weix',
    -DBNAME => "ensembl_website_$version",
    -PASS => 'warelab',
    -SPECIES => 'multi',
    -GROUP => 'web'
  );
  Bio::EnsEMBL::Production::DBSQL::DBAdaptor->new(
    -HOST => 'cabot',
    -PORT => 3306,
    -USER => 'weix',
    -DBNAME => 'ensembl_production',
    -PASS => 'warelab',
    -SPECIES => 'multi',
    -GROUP => 'production'
  );
}
1;
