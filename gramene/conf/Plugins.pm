## If you wish to use the EnsEMBL web-code from the command line, you will
## need to hardcode the server root here 

## $SiteDefs::ENSEMBL_SERVERROOT = '/path to root of ensembl tree';
my $ENSROOT = $SiteDefs::ENSEMBL_SERVERROOT;
my $GRMENSROOT = $SiteDefs::ENSEMBL_SERVERROOT.'/gramene-live';
$SiteDefs::ENSEMBL_PLUGINS = [
  'EnsEMBL::Gramene'       => $GRMENSROOT/gramene',
  'EnsEMBL::Mart'          => $ENSROOT.'/public-plugins/mart',
  'EG::Plants' => $ENSROOT.'/eg-plugins/plants',
  'EG::Common'=> $ENSROOT.'/eg-plugins/common',
#  'EnsEMBL::Memcached' => $ENSROOT.'/public-plugins/memcached',
];

1;
