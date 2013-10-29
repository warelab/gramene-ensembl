## If you wish to use the EnsEMBL web-code from the command line, you will
## need to hardcode the server root here 

## $SiteDefs::ENSEMBL_SERVERROOT = '/path to root of ensembl tree';
## Copy this file to the $ENSROOT/conf file,
## Remove $ENSROOT/conf/config.packed, and
## Restart the server

my $ENSROOT = $SiteDefs::ENSEMBL_SERVERROOT;
my $GRMROOT = $SiteDefs::ENSEMBL_SERVERROOT.'/gramene-live';
$SiteDefs::ENSEMBL_PLUGINS = [
  'EnsEMBL::MaizeBlast'  => $GRMROOT.'/ensembl-plugins/maizeblast',
  'EnsEMBL::Maize'       => $GRMROOT.'/ensembl-plugins/maize',
  'EnsEMBL::GrmBlast'    => $GRMROOT.'/ensembl-plugins/grmblast',
#  'EnsEMBL::Gramene'     => $GRMROOT.'/ensembl-plugins/gramene',
];

1;
