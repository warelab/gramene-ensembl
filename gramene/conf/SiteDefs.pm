package EnsEMBL::Gramene::SiteDefs;
use strict;

# These are the gramene-specific edits to the main Ensembl SiteDefs.pm file
sub update_conf {

   #$SiteDefs::ENSEMBL_SERVERNAME             = 'ensembl.gramene.org';
  $SiteDefs::ENSEMBL_SERVERNAME             = 'ensembl-dev.gramene.org';

  $SiteDefs::ENSEMBL_MAX_PROCESS_SIZE     = 2000000; 
  $SiteDefs::ENSEMBL_BASE_URL     = $SiteDefs::ENSEMBL_SERVERNAME;
  $SiteDefs::SITE_RELEASE_VERSION = 68;
  $SiteDefs::SITE_RELEASE_VERSION_EG = 57;
  $SiteDefs::SITE_RELEASE_DATE    = 'Nov 2023';  #G66 is 'Dec 2022'
  $SiteDefs::SITE_NAME            = 'Gramene';
  $SiteDefs::ENSEMBL_SITETYPE = 'Ensembl Plants';
  $SiteDefs::SITE_FTP             = 'https://ftp.gramene.org/pub';  #the one that's used for download link on species page is ENSEMBL_FTP_URL in DEFAULT.ini
  $SiteDefs::PE_URL             = 'https://plants.ensembl.org';
  $SiteDefs::ENSEMBL_PORT       = 80;
  $SiteDefs::ENSEMBL_PROXY_PORT = 80; # Port used for self-referential URLs
  $SiteDefs::ENSEMBL_USER       = 'nobody';#getpwuid($>);          
  $SiteDefs::ENSEMBL_GROUP      = 'nobody';#getgrgid($));           

  $SiteDefs::ENSEMBL_SERVERADMIN            = 'weix@cshl.edu';
  $SiteDefs::ENSEMBL_MAIL_SERVER            = 'localhost';

  $SiteDefs::MAX_PROCESS_SIZE = 2000000;
  $SiteDefs::GRAPHIC_TTF_PATH     = "/usr/share/fonts/msttcorefonts/";


  $SiteDefs::SAMTOOLS_DIR = $SiteDefs::ENSEMBL_SERVERROOT.'/samtools'; 

  $SiteDefs::ENSEMBL_DEBUG_FLAGS             = 0; # 24;
  $SiteDefs::ENSEMBL_LONGPROCESS_MINTIME     = 10;

  $SiteDefs::DATAFILE_ROOT        = '/usr/local';                                  ## Base path for ro data files
  $SiteDefs::DATAFILE_BASE_PATH   = '/usr/local/vcf';


  $SiteDefs::EBEYE_REST_ENDPOINT  = 'https://data.gramene.org/ebeye' . $SiteDefs::SITE_RELEASE_VERSION;
   #'https://data.gramene.org/ebeye66';

  #$SiteDefs::NCBIBLAST_REST_ENDPOINT = 'http://brie:5202';
###############################################################################
# this section copied over from ensembl-webcode/conf/SiteDef.pm
## GDPR variables
## Some variables are assigned null for external users to override
###############################################################################
$SiteDefs::GDPR_VERSION                 = 1;
$SiteDefs::GDPR_COOKIE_NAME             = 'gdpr';
#our $GDPR_POLICY_URL              = 'https://www.ebi.ac.uk/data-protection/ensembl/privacy-notice';
#our $GDPR_TERMS_URL               = 'https://www.ebi.ac.uk/about/terms-of-use';


  $SiteDefs::ENSEMBL_TOOLS_PERL_BIN         = '/usr/local/bin/perl';
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
  $SiteDefs::ENSEMBL_CHAIN_FILE_DIR       = '/usr/local/ensembl-live/tools_data/assembly_chain';

$SiteDefs::ENSEMBL_VEP_CACHE_DIR              = "/usr/local/ensembl-live/tools_data/vep/";
$SiteDefs::ENSEMBL_VEP_PLUGIN_DATA_DIR        = "/usr/local/ensembl-live/tools_data/vep/Plugins";                       # path to vep plugin data files on the LSF host (or local machine if job running locally) 
 
push @$SiteDefs::ENSEMBL_EXTRA_INC, '/usr/local/ensembl-live/htslib', '/usr/local/BioPerl-1.6.922', '/usr/local/ensembl-live/ensembl-io/modules', '/usr/local/ensembl-live/ensembl-funcgen/modules', '/usr/local/ensembl-live/ensembl-variation/modules';

 
$SiteDefs::ENSEMBL_VEP_SCRIPT_DEFAULT_OPTIONS = {                                                 # Default options for command line vep script (keys with value undef get ignored)
    'host'        => 'colden',                                                                       # Database host (defaults to ensembldb.ensembl.org)
    'user'        => 'weix',                                                                       # Defaults to 'anonymous'
    'password'    => 'warelab',                                                                       # Not used by default
    'port'        => 3306,                                                                       # Defaults to 5306
    'fork'        => 4,                                                                           # Enable forking, using 4 forks
  };

  #----------
  # User database
  $SiteDefs::ENSEMBL_USERDB_NAME = 'ensembl_accounts';  #changed to ensembl_accounts lately
  $SiteDefs::ENSEMBL_USERDB_USER = 'gramene_web';
  $SiteDefs::ENSEMBL_USERDB_HOST = 'colden.cshl.edu';
  $SiteDefs::ENSEMBL_USERDB_PORT =  3306;
  $SiteDefs::ENSEMBL_USERDB_PASS = 'gram3n3';

  #----------
  # Logging
  $SiteDefs::ENSEMBL_LOGDIR = $SiteDefs::ENSEMBL_SERVERROOT."/logs";
  $SiteDefs::ENSEMBL_PIDFILE   = "$SiteDefs::ENSEMBL_LOGDIR/httpd.pid";
  $SiteDefs::ENSEMBL_ERRORLOG  = "$SiteDefs::ENSEMBL_LOGDIR/error.log";
  $SiteDefs::ENSEMBL_CUSTOMLOG = "$SiteDefs::ENSEMBL_LOGDIR/access.log combined";
#our $ENSEMBL_LOGDIR               = defer { $SiteDefs::ENSEMBL_SERVERROOT."/logs" };              # Path for log files, used to be $ENSEMBL_SYS_DIR/logs/$ENSEMBL_SERVER_SIGNATURE
#our $ENSEMBL_PIDFILE              = defer { "$ENSEMBL_LOGDIR/httpd.pid" };                                    # httpd process id
#our $ENSEMBL_ERRORLOG             = defer { "$ENSEMBL_LOGDIR/error_log" };                                    # Error log file
#our $ENSEMBL_CUSTOMLOG            = defer { "$ENSEMBL_LOGDIR/access_log ensembl_extended" };    


  #----------
  # Mart/Blast
  $SiteDefs::ENSEMBL_BLAST_ENABLED = 1; # Creates header link for blast
  $SiteDefs::ENSEMBL_BLAST_BY_SEQID = 1; # blast on gene page sequence
  $SiteDefs::ENSEMBL_MART_ENABLED = 1; # And mart

  $SiteDefs::ENSEMBL_VEP_ENABLED    = 1;
  $SiteDefs::ENSEMBL_AC_ENABLED     = 1;
  $SiteDefs::ENSEMBL_IDM_ENABLED    = 1;
  $SiteDefs::ENSEMBL_FC_ENABLED     = 0;
  $SiteDefs::ENSEMBL_LD_ENABLED     = 0;
  $SiteDefs::ENSEMBL_VP_ENABLED     = 0;
  $SiteDefs::ENSEMBL_DS_ENABLED     = 0;

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
    # 'Zea_mays'             => [('zea_mays', 'zm','maize')],
    # 'Physcomitrella_patens' => [('physcomitrella_patens', 'pp','physcomitrella')],
    # 'Arabidopsis_thaliana'  => [('arabidopsis_thaliana')],
    # );
#  push(@{$SiteDefs::__species_aliases{'Populus_trichocarpa'}},'poplar');

$SiteDefs::ENSEMBL_TOOLS_LIST = [
    'Blast'             => 'BLAST/BLAT',
    'VEP'               => 'Variant Effect Predictor',
    #'FileChameleon'     => 'File Chameleon',
    'AssemblyConverter' => 'Assembly Converter',
    'IDMapper'          => 'ID History Converter',
    #'AlleleFrequency'   => 'Allele Frequency Calculator',
    #'VcftoPed'          => 'VCF to PED Converter',
    #'DataSlicer'        => 'Data Slicer',
    #'VariationPattern'  => 'Variation Pattern Finder',
    #'LD'                => 'Linkage Disequilibrium Calculator',
  ];


}

1;
