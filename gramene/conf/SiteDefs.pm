package EnsEMBL::Gramene::SiteDefs;
use strict;

# These are the gramene-specific edits to the main Ensembl SiteDefs.pm file
sub update_conf {

  $SiteDefs::ENSEMBL_SERVERNAME             = 'ensembl.gramene.org';

  $SiteDefs::ENSEMBL_BASE_URL     = $SiteDefs::ENSEMBL_SERVERNAME;
  $SiteDefs::SITE_RELEASE_VERSION = "42"; 
  $SiteDefs::SITE_RELEASE_VERSION_EG = "23";
  #$SiteDefs::ENSEMBL_VERSION = 65;
  $SiteDefs::SITE_RELEASE_DATE    = 'May 2014';
  $SiteDefs::SITE_NAME            = 'Gramene';
  $SiteDefs::SITE_FTP             = 'ftp://ftp.gramene.org/pub';
  $SiteDefs::GRAMENE_FTP_URL	  = 'ftp://ftp.gramene.org/pub';
  $SiteDefs::PE_URL             = 'http://plants.ensembl.org';
  $SiteDefs::ENSEMBL_PORT       = 80;
  $SiteDefs::ENSEMBL_PROXY_PORT = 80; # Port used for self-referential URLs
  $SiteDefs::ENSEMBL_USER       = 'nobody';#getpwuid($>);          
  $SiteDefs::ENSEMBL_GROUP      = 'nobody';#getgrgid($));           

  $SiteDefs::SAMTOOLS_DIR = $SiteDefs::ENSEMBL_SERVERROOT.'/samtools'; 

  #----------
  # User database
  $SiteDefs::ENSEMBL_USERDB_NAME = 'ensembl_web_user_db';
  #'ensembl_web_user_db_31_57';
  $SiteDefs::ENSEMBL_USERDB_USER = 'gramene_web';
  $SiteDefs::ENSEMBL_USERDB_HOST = 'colden.cshl.edu';
  $SiteDefs::ENSEMBL_USERDB_PORT =  3306;
  $SiteDefs::ENSEMBL_USERDB_PASS = 'gram3n3';

  #----------
  # Logging
  #my $LOG_ROOT = $SiteDefs::ENSEMBL_WEBROOT."/logs";
  #$SiteDefs::ENSEMBL_PIDFILE   = "$LOG_ROOT/httpd.pid";
  #$SiteDefs::ENSEMBL_ERRORLOG  = "$LOG_ROOT/error.log";
  #$SiteDefs::ENSEMBL_CUSTOMLOG = "$LOG_ROOT/access.log combined";

  #----------
  # Mart/Blast
  $SiteDefs::ENSEMBL_BLAST_ENABLED ++; # Creates header link for blast
  $SiteDefs::ENSEMBL_MART_ENABLED ++; # And mart
  push @SiteDefs::ENSEMBL_HTDOCS_DIRS, # Needed due to EG plugin
    $SiteDefs::ENSEMBL_SERVERROOT.'/biomart-perl/htdocs';

  #----------
  # Temp files
  #$SiteDefs::ENSEMBL_TMP_URL = "/ens-tmp";

  #----------
  #GeoLiteCity database file
  $SiteDefs::GEOCITY_DAT = $SiteDefs::ENSEMBL_SERVERROOT.'/geocity/GeoLiteCity.dat';

  #----------
  # Paths and so on
  # This should be set by the plugin, but not working. Need to hard code!
#  push @SiteDefs::ENSEMBL_PERL_DIRS, $SiteDefs::ENSEMBL_SERVERROOT 
#      .'/gramene-live/ensembl-plugins/gramene/perl';

  #----------
  # Species stuff
  #push @$SiteDefs::ENSEMBL_DATASETS, 'Oryza_barthii';
  #push @$SiteDefs::ENSEMBL_DATASETS, 'Oryza_nivara';
  push @$SiteDefs::ENSEMBL_DATASETS, 'Oryza_rufipogon';
  push @$SiteDefs::ENSEMBL_DATASETS, 'Oryza_rufipogon3s';
  push @$SiteDefs::ENSEMBL_DATASETS, 'Oryza_minutabb';
  push @$SiteDefs::ENSEMBL_DATASETS, 'Oryza_minutacc';
  push @$SiteDefs::ENSEMBL_DATASETS, 'Oryza_officinalis';
  #push @$SiteDefs::ENSEMBL_DATASETS, 'Oryza_punctata';
  #push @$SiteDefs::ENSEMBL_DATASETS, 'Oryza_glumaepatula';
  #push @$SiteDefs::ENSEMBL_DATASETS, 'Oryza_meridionalis'; 
  push @$SiteDefs::ENSEMBL_DATASETS, 'Oryza_granulata';
  push @$SiteDefs::ENSEMBL_DATASETS, 'Oryza_longistaminata';
  push @$SiteDefs::ENSEMBL_DATASETS, 'Leersia_perrieri';

#warn("species datases are", @$SiteDefs::ENSEMBL_DATASETS);
 #push @$SiteDefs::ENSEMBL_DATASETS, 'Zea_mays';
  #push @$SiteDefs::ENSEMBL_DATASETS, 'Physcomitrella_patens';

  $SiteDefs::ENSEMBL_PRIMARY_SPECIES    = 'Oryza_sativa'; # Default
  $SiteDefs::ENSEMBL_SECONDARY_SPECIES  = 'Arabidopsis_thaliana';
  %SiteDefs::__species_aliases =
    (
     # These are supplimental species to EnsemblGenomes
     %SiteDefs::__species_aliases,

     'Leersia_perrieri'    => [('lper','leersia')],
   #  'Oryza_barthii'        => [('ob','barthii')],
     'Oryza_brachyantha'    => [('obr','brachyantha')],
     'Oryza_glaberrima3s'     => [('og3s','glaberrima3s')],
#     'Oryza_glaberrima'     => [('og','glaberrima')],
    # 'Oryza_glumaepatula'    => [('oglu','glumaepatula')],
     'Oryza_granulata'    => [('ogra','granulata')],
     'Oryza_indica'         => [('oryza_sativa_indica')],
     'Oryza_longistaminata'    => [('olon','longistaminata')],
    # 'Oryza_meridionalis'       => [('omer','meridionalis')],
     'Oryza_minutabb'       => [('omBB','Oryza_minutabb', 'Oryza_minutaBB','minutaBB')],
     'Oryza_minutacc'       => [('omCC', 'Oryza_minutacc', 'Oryza_minutaCC','minutaCC')],
    # 'Oryza_nivara'         => [('on','nivara')],
     'Oryza_officinalis'    => [('oo','officianlis')],
    # 'Oryza_punctata'       => [('op','punctata')],
     'Oryza_rufipogon3s'    => [('or3s', 'oryza_rufipogon')],
     'Oryza_rufipogon'      => [('or','rufipogon','oryza_rufipogon_fpc')],
     'Oryza_sativa'         => [('oryza_sativa_japonica')],

     # Needed for the home page and species lists
     #'Zea_mays'             => [('zm','maize')],
     #'Physcomitrella_patens' => [('pp','physcomitrella')],
     );

#  push(@{$SiteDefs::__species_aliases{'Oryza_indica'}},'Oryza_sativa_indica');
#  push(@{$SiteDefs::__species_aliases{'Oryza_indica'}},'indica');
#  push(@{$SiteDefs::__species_aliases{'Oryza_sativa'}},'Oryza_sativa_japonica');
#  push(@{$SiteDefs::__species_aliases{'Populus_trichocarpa'}},'poplar');
}

1;
