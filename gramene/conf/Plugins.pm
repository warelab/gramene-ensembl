## If you wish to use the EnsEMBL web-code from the command line, you will
## need to hardcode the server root here 

## $SiteDefs::ENSEMBL_SERVERROOT = '/path to root of ensembl tree';
my $ENSROOT = $SiteDefs::ENSEMBL_SERVERROOT;
my $GRMROOT = $SiteDefs::ENSEMBL_SERVERROOT.'/gramene-live';
$SiteDefs::ENSEMBL_PLUGINS = [
  'EnsEMBL::Gramene'       => $GRMROOT.'/gramene',
  'EnsEMBL::Mart'          => $ENSROOT.'/public-plugins/mart',
  'EG::EBEyeSearch'        => $ENSROOT.'/eg-web-search/',
  'EnsEMBL::Genoverse'     => $ENSROOT.'/public-plugins/genoverse',
  'EG::Plants'             => $ENSROOT.'/eg-web-plants/',
  'EG::Common'             => $ENSROOT.'/eg-web-common/',
  'EG::Widgets' 	   => $ENSROOT.'/public-plugins/widgets',
  'EnsEMBL::Tools'         => $ENSROOT.'/public-plugins/tools',
  'EnsEMBL::Users'         => $ENSROOT.'/public-plugins/users',
  'EnsEMBL::Memcached'     => $ENSROOT.'/public-plugins/memcached',
  'EnsEMBL::Doc'           => $ENSROOT.'/public-plugins/docs',
];

1;
