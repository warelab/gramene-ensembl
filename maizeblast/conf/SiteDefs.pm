package EnsEMBL::MaizeBlast::SiteDefs;
use strict;

# These are the gramene-specific edits to the main Ensembl SiteDefs.pm file
sub update_conf {
  warn( "===> ", __PACKAGE__, "->update_conf\n" );
  my $DEVELOPMENT_SITE       = 1;    #acorn or filetta

  $SiteDefs::ENSEMBL_PORT         = 8028;
  $SiteDefs::ENSEMBL_PROXY_PORT   = 80; # Port used for self-referential URLs:
  $SiteDefs::ENSEMBL_USER         = 'nobody';     
  $SiteDefs::ENSEMBL_GROUP         = 'nobody';
  $SiteDefs::ENSEMBL_SERVERNAME   = 'www.maizesequence.org';
  $SiteDefs::ENSEMBL_MAIL_ERRORS  = 1;
  $SiteDefs::ENSEMBL_ERRORS_TO    = 'webmaster@maizesequence.org';

  #----------
  # Ensembl user_db connection info
  $SiteDefs::ENSEMBL_USERDB_NAME  = 'ensembl_web_user_db';
  $SiteDefs::ENSEMBL_USERDB_USER  = 'maize_rw';
  $SiteDefs::ENSEMBL_USERDB_HOST  = 'cannon.cshl.edu';
  $SiteDefs::ENSEMBL_USERDB_PORT  = '3306';
  $SiteDefs::ENSEMBL_USERDB_PASS  = 'z3@m@ys';


  #----------
  # tmp dirs; tweak for proxy
  $SiteDefs::ENSEMBL_TMP_URL       = '/BlastView/tmp';
  $SiteDefs::ENSEMBL_TMP_URL_IMG   = '/BlastView/img-tmp';
  $SiteDefs::ENSEMBL_TMP_URL_CACHE = '/BlastView/img-cache';

  $SiteDefs::ENSEMBL_PRIMARY_SPECIES   = 'Zea_mays2';
  $SiteDefs::ENSEMBL_SECONDARY_SPECIES = 'Zea_mays2';
  %SiteDefs::__species_aliases = (
        # 'Zea_mays'  => [ ('maize',        'zm') ],
	'Zea_mays2' => [ ('zea_mays_bac', 'maize_bac', 'zm2' ) ],
  );
}

1;
