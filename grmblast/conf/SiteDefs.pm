package EnsEMBL::GrmBlast::SiteDefs;
use strict;

# These are the Gramene BlastView-specific edits to the main Gramene
# SiteDefs.pm file
sub update_conf {
  warn( "===> ", __PACKAGE__, "->update_conf\n" );

  #$SiteDefs::ENSEMBL_PORT       = '8018';
  $SiteDefs::ENSEMBL_PORT       = '8019';

  #----------
  # tmp dirs; tweak for proxy
  $SiteDefs::ENSEMBL_TMP_URL       = '/BlastView/tmp';
  $SiteDefs::ENSEMBL_TMP_URL_IMG   = '/BlastView/img-tmp';
  $SiteDefs::ENSEMBL_TMP_URL_CACHE = '/BlastView/img-cache';
  $SiteDefs::ENSEMBL_TMP_DIR_MINIFIED = '/BlastView/minified'; 
      #= $SiteDefs::ENSEMBL_WEBROOT.'/htdocs/BlastView/minified';
  $SiteDefs::ENSEMBL_TMP_URL_MINIFIED = '/BlastView/minified';

  #----------
  # Species stuff
  # Ditch problematic species
  #delete( $SiteDefs::__species_aliases{Oryza_rufipogon} );
  #delete( $SiteDefs::__species_aliases{Zea_mays} );
  #delete( $SiteDefs::__species_aliases{Zea_mays2} );
}

1;
