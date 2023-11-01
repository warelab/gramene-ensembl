## If you wish to use the EnsEMBL web-code from the command line, you will
## need to hardcode the server root here 

## $SiteDefs::ENSEMBL_SERVERROOT = '/path to root of ensembl tree';
## Copy this file to the $ENSROOT/conf file,
## Remove $ENSROOT/conf/config.packed, and
## Restart the server

my $ENSROOT = $SiteDefs::ENSEMBL_SERVERROOT;
my $GRMROOT = $SiteDefs::ENSEMBL_SERVERROOT.'/gramene-live';
$SiteDefs::ENSEMBL_PLUGINS = [
  # grmblast contains the ini files to configure BlastView
  'EnsEMBL::GrmBlast'  => $GRMROOT.'/grmblast',
  'EnsEMBL::Gramene'   => $GRMROOT.'/gramene',
  'EG::Plants' => $ENSROOT.'/eg-web-plants/',
  'EG::Common'=> $ENSROOT.'/eg-web-common/',
];

1;