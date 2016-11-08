package EnsEMBL::Gramene::SiteDefs;
use strict;

# These are the gramene-specific edits to the main Ensembl SiteDefs.pm file
sub update_conf {

  $SiteDefs::ENSEMBL_SERVERNAME             = 'ensembl.gramene.org';
#  $SiteDefs::ENSEMBL_SERVERNAME             = 'dev.gramene.org';

  $SiteDefs::ENSEMBL_BASE_URL     = $SiteDefs::ENSEMBL_SERVERNAME;
  $SiteDefs::SITE_RELEASE_VERSION = 52; 
  $SiteDefs::SITE_RELEASE_VERSION_EG = 33;
  #$SiteDefs::ENSEMBL_VERSION = 65;
  $SiteDefs::SITE_RELEASE_DATE    = 'Nov 2016';
  $SiteDefs::SITE_NAME            = 'Gramene';
  $SiteDefs::ENSEMBL_SITETYPE = 'Gramene';
  $SiteDefs::SITE_FTP             = 'ftp://ftp.gramene.org/pub';
  $SiteDefs::GRAMENE_FTP_URL	  = 'ftp://ftp.gramene.org/pub';
#  $SiteDefs::OGE_FTP_URL          = 'ftp://ftp.gramene.org/pub/gramene/oge/release-current'
  $SiteDefs::PE_URL             = 'http://plants.ensembl.org';
  $SiteDefs::ENSEMBL_PORT       = 80;
  $SiteDefs::ENSEMBL_PROXY_PORT = 80; # Port used for self-referential URLs
  $SiteDefs::ENSEMBL_USER       = 'nobody';#getpwuid($>);          
  $SiteDefs::ENSEMBL_GROUP      = 'nobody';#getgrgid($));           

  $SiteDefs::ENSEMBL_SERVERADMIN            = 'weix@cshl.edu';
  $SiteDefs::ENSEMBL_MAIL_SERVER            = 'localhost';

  $SiteDefs::SAMTOOLS_DIR = $SiteDefs::ENSEMBL_SERVERROOT.'/samtools'; 

  $SiteDefs::ENSEMBL_DEBUG_FLAGS             = 0; # 24;
  $SiteDefs::ENSEMBL_LONGPROCESS_MINTIME     = 10;

  $SiteDefs::EBEYE_REST_ENDPOINT     = 'http://data.gramene.org/ebeye';

  $SiteDefs::ENSEMBL_TOOLS_PERL_BIN             = '/usr/local/bin/perl';
  $SiteDefs::ENSEMBL_TMP_DIR_BLAST          = $SiteDefs::ENSEMBL_SERVERROOT."/blastqueue";
  $SiteDefs::ENSEMBL_BLASTSCRIPT            = $SiteDefs::ENSEMBL_WEBROOT."/utils/runblast.pl";

#  $SiteDefs::ENSEMBL_BLAT_BIN_PATH          = '/usr/local/ucscblat/blat';
#  $SiteDefs::ENSEMBL_BLAT_TWOBIT_DIR        = '/usr/local/ucscblat/';
  $SiteDefs::ENSEMBL_NCBIBLAST_BIN_PATH     = '/usr/local/ncbi-blast-2.2.30+/bin'; # path to blast executables  
#  $SiteDefs::ENSEMBL_NCBIBLAST_MATRIX       = '/path/to/ncbi-blast/data'; # path to blast matrix files 
  $SiteDefs::ENSEMBL_NCBIBLAST_DATA_PATH_DNA = "/usr/local/blastdb/ncbi_blast/genomic"; # path for the blast DNA index files 
  $SiteDefs::ENSEMBL_NCBIBLAST_DATA_PATH    = "/usr/local/blastdb/ncbi_blast/genes"; # path for the blast index files (other than DNA) 
  $SiteDefs::ENSEMBL_REPEATMASK_BIN_PATH    = '/usr/local/RepeatMasker'; # path to RepeatMasker executable
  $SiteDefs::ASSEMBLY_CONVERTER_BIN_PATH = '/usr/local/bin/CrossMap.py';
  $SiteDefs::ENSEMBL_CHAIN_FILE_DIR       = '/usr/local/ensembl-live/tools_data/assembly_converter';
  #$SiteDefs::IDMAPPER_SCRIPT              = '/usr/local/ensembl-live/ensembl-tools/scripts/id_history_converter/IDmapper.pl'; 
 
$SiteDefs::ENSEMBL_VEP_CACHE_DIR              = "/usr/local/ensembl-live/tools_data/vep/";
$SiteDefs::ENSEMBL_VEP_PLUGIN_DATA_DIR        = "/usr/local/ensembl-live/tools_data/vep/Plugins";                       # path to vep plugin data files on the LSF host (or local machine if job running locally) 
  
$SiteDefs::ENSEMBL_VEP_SCRIPT_DEFAULT_OPTIONS = {                                                 # Default options for command line vep script (keys with value undef get ignored)
    '--host'        => 'colden',                                                                       # Database host (defaults to ensembldb.ensembl.org)
    '--user'        => 'weix',                                                                       # Defaults to 'anonymous'
    '--password'    => 'warelab',                                                                       # Not used by default
    '--port'        => 3306,                                                                       # Defaults to 5306
    '--fork'        => 4,                                                                           # Enable forking, using 4 forks
  };


  #----------
  # User database
  $SiteDefs::ENSEMBL_USERDB_NAME = 'ensembl_accounts_86';  #changed to ensembl_accounts lately
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
  $SiteDefs::ENSEMBL_BLAST_ENABLED = 1; # Creates header link for blast
  $SiteDefs::ENSEMBL_MART_ENABLED = 1; # And mart

  $SiteDefs::ENSEMBL_VEP_ENABLED    = 1;
  $SiteDefs::ENSEMBL_AC_ENABLED     = 1;
  $SiteDefs::ENSEMBL_IDM_ENABLED    = 0;

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

  $SiteDefs::ENSEMBL_PRIMARY_SPECIES    = 'Oryza_sativa'; # Default
  $SiteDefs::ENSEMBL_SECONDARY_SPECIES  = 'Arabidopsis_thaliana';
  #%SiteDefs::__species_aliases =
   # (
     # These are supplimental species to EnsemblGenomes
    # %SiteDefs::__species_aliases,
     #'Zea_mays'             => [('zm','maize')],
     #'Physcomitrella_patens' => [('pp','physcomitrella')],
     #);
#  push(@{$SiteDefs::__species_aliases{'Populus_trichocarpa'}},'poplar');

}

1;
