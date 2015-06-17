package EnsEMBL::Gramene::SiteDefs;
use strict;

# These are the gramene-specific edits to the main Ensembl SiteDefs.pm file
sub update_conf {

  $SiteDefs::ENSEMBL_SERVERNAME             = 'oge.gramene.org';

  $SiteDefs::ENSEMBL_BASE_URL     = $SiteDefs::ENSEMBL_SERVERNAME;
  $SiteDefs::SITE_RELEASE_VERSION = 2; 
  $SiteDefs::SITE_RELEASE_DATE    = 'Jun 2015';
  $SiteDefs::SITE_NAME            = 'Gramene OGE';
  $SiteDefs::ENSEMBL_SITETYPE     = 'Gramene OGE';
  $SiteDefs::SITE_FTP             = 'ftp://ftp.gramene.org/pub';
  $SiteDefs::GRAMENE_FTP_URL	  = 'ftp://ftp.gramene.org/pub';
#  $SiteDefs::OGE_FTP_URL          = 'ftp://ftp.gramene.org/pub/gramene/oge/release-current';
  $SiteDefs::PE_URL             = 'http://plants.ensembl.org';
  $SiteDefs::ENSEMBL_PORT       = 80;
  $SiteDefs::ENSEMBL_PROXY_PORT = 80; # Port used for self-referential URLs
  $SiteDefs::ENSEMBL_USER       = 'nobody';#getpwuid($>);          
  $SiteDefs::ENSEMBL_GROUP      = 'nobody';#getgrgid($));           

  $SiteDefs::ENSEMBL_SERVERADMIN            = 'weix@cshl.edu';
  $SiteDefs::ENSEMBL_MAIL_SERVER       = 'localhost';
  $SiteDefs::ENSEMBL_SERVERNAME        = 'oge.gramene.org';
  

  $SiteDefs::SAMTOOLS_DIR = $SiteDefs::ENSEMBL_SERVERROOT.'/samtools'; 

  $SiteDefs::ENSEMBL_DEBUG_FLAGS             = 24;
  $SiteDefs::ENSEMBL_LONGPROCESS_MINTIME     = 10;

  $SiteDefs::ENSEMBL_TMP_DIR_BLAST          = $SiteDefs::ENSEMBL_SERVERROOT."/blastqueue";
  $SiteDefs::ENSEMBL_BLASTSCRIPT            = $SiteDefs::ENSEMBL_WEBROOT."/utils/runblast.pl";

# if this is a blast tool web site you need to set the following
# All of the Ensembl tools (blast tool here) use some kind of index/cache files for fast data processing, 
# and it makes sense to keep all of them in the same place for easy maintenance. 
# Create this directory somewhere that the webserver can access

   $SiteDefs::ENSEMBL_NCBIBLAST_BIN_PATH     = '/usr/local/ncbi-blast-2.2.30+/bin'; # path to blast executables  
   $SiteDefs::ENSEMBL_NCBIBLAST_DATA_PATH_DNA = "/usr/local/blastdb/ncbi_blast/genomic"; # path for the blast DNA index files 
   $SiteDefs::ENSEMBL_NCBIBLAST_DATA_PATH    = "/usr/local/blastdb/ncbi_blast/genes"; # path for the blast index files (other than DNA) 
   $SiteDefs::ENSEMBL_REPEATMASK_BIN_PATH    = '/usr/local/RepeatMasker'; # path to RepeatMasker executable


  #----------
  # User database
  $SiteDefs::ENSEMBL_USERDB_NAME = 'ensembl_accounts';
  #'ensembl_web_user_db_31_57';
  $SiteDefs::ENSEMBL_USERDB_USER = 'gramene_web';
  $SiteDefs::ENSEMBL_USERDB_HOST = 'cabot.cshl.edu';
  $SiteDefs::ENSEMBL_USERDB_PORT =  3306;
  $SiteDefs::ENSEMBL_USERDB_PASS = 'gram3n3';

  #----------
  # Mart/Blast
  $SiteDefs::ENSEMBL_MART_ENABLED   = 0; # And mart, do we need mart for OGE?
  $SiteDefs::ENSEMBL_BLAST_ENABLED  = 1;
  $SiteDefs::ENSEMBL_VEP_ENABLED    = 0;
  $SiteDefs::ENSEMBL_AC_ENABLED     = 0;

  $SiteDefs::ENSEMBL_TOOLS_JOB_DISPATCHER = { 'Blast' => 'NcbiBlast' };
  $SiteDefs::NCBIBLAST_REST_ENDPOINT  = 'http://www.ebi.ac.uk/Tools/services/rest/ncbiblast';


  push @SiteDefs::ENSEMBL_HTDOCS_DIRS, # Needed due to EG plugin
    $SiteDefs::ENSEMBL_SERVERROOT.'/biomart-perl/htdocs';


  #----------
  #GeoLiteCity database file
  $SiteDefs::GEOCITY_DAT = $SiteDefs::ENSEMBL_SERVERROOT.'/geocity/GeoLiteCity.dat';


  #----------
  # Species stuff

  $SiteDefs::ENSEMBL_DATASETS = undef;	
  $SiteDefs::ENSEMBL_DATASETS = [sort qw(
        Arabidopsis_thaliana
        Brachypodium_distachyon
        Leersia_perrieri
        Oryza_barthii
        Oryza_brachyantha
        Oryza_glaberrima
        Oryza_glumaepatula
        Oryza_granulata3s
        Oryza_indica
        Oryza_longistaminata
        Oryza_meridionalis
        Oryza_minutabb3s
        Oryza_minutacc3s
        Oryza_nivara
        Oryza_officinalis3s
        Oryza_punctata
        Oryza_rufipogon
        Oryza_sativa
        Sorghum_bicolor
    )];

  $SiteDefs::ENSEMBL_PRIMARY_SPECIES    = 'Oryza_sativa'; # Default
  $SiteDefs::ENSEMBL_SECONDARY_SPECIES  = 'Arabidopsis_thaliana';

}

1;
