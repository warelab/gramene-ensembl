package EnsEMBL::GrameneDev::SiteDefs;
use strict;

# These are the gramene development site modifications to the live site
# config. Mainly hostname differences

sub update_conf {
  #warn( "===> ", __PACKAGE__, "->update_conf\n" );
  $SiteDefs::ENSEMBL_SERVERNAME   = 'dev.gramene.org';  
  $SiteDefs::ENSEMBL_USERDB_HOST  = 'cabot.cshl.edu';
}
1;
