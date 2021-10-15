package EnsEMBL::Gramene::SiteDefs;
use strict;

# These are the gramene-specific edits to the main Ensembl SiteDefs.pm file
sub update_conf {

  $SiteDefs::ENSEMBL_MAX_PROCESS_SIZE     = 2 * 1024 * 1024;
  $SiteDefs::ENSEMBL_SERVERNAME             = 'oge-ensembl.gramene.org';

  $SiteDefs::ENSEMBL_BASE_URL     = $SiteDefs::ENSEMBL_SERVERNAME;
  $SiteDefs::SITE_RELEASE_VERSION = 3; 
  $SiteDefs::SITE_RELEASE_DATE    = 'Aug 2018';
  $SiteDefs::SITE_NAME            = 'Gramene OGE';
  $SiteDefs::ENSEMBL_SITETYPE     = 'Gramene OGE';
  $SiteDefs::SITE_FTP             = 'ftp://ftp.gramene.org/pub';
  $SiteDefs::GRAMENE_FTP_URL	  = 'ftp://ftp.gramene.org/pub';
#  $SiteDefs::OGE_FTP_URL          = 'ftp://ftp.gramene.org/pub/gramene/oge/release-current';
  $SiteDefs::PE_URL             = 'http://plants.ensembl.org';
  $SiteDefs::ENSEMBL_PORT       = 86;
  $SiteDefs::ENSEMBL_PROXY_PORT = 86; # Port used for self-referential URLs
  $SiteDefs::ENSEMBL_USER       = 'nobody';#getpwuid($>);          
  $SiteDefs::ENSEMBL_GROUP      = 'nobody';#getgrgid($));           

  $SiteDefs::ENSEMBL_SERVERADMIN            = 'weix@cshl.edu';
  $SiteDefs::ENSEMBL_MAIL_SERVER       = 'localhost';

  $SiteDefs::SAMTOOLS_DIR = $SiteDefs::ENSEMBL_SERVERROOT.'/samtools'; 

  $SiteDefs::ENSEMBL_DEBUG_FLAGS             = 0;
  $SiteDefs::ENSEMBL_LONGPROCESS_MINTIME     = 10;

$SiteDefs::EBEYE_REST_ENDPOINT     = 'http://data.gramene.org/ebeye';

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
  #$SiteDefs::ENSEMBL_USERDB_NAME = 'ensembl_accounts_oge3';
  #'ensembl_web_user_db_31_57';
  $SiteDefs::ENSEMBL_USERDB_USER = 'gramene_web';
  $SiteDefs::ENSEMBL_USERDB_HOST = 'cabot.cshl.edu';
  $SiteDefs::ENSEMBL_USERDB_PORT =  3306;
  $SiteDefs::ENSEMBL_USERDB_PASS = 'gram3n3';

  #----------
  # Mart/Blast
  $SiteDefs::ENSEMBL_MART_ENABLED   = 0; # And mart, do we need mart for OGE?
  $SiteDefs::ENSEMBL_BLAST_ENABLED  = 0;
  $SiteDefs::ENSEMBL_VEP_ENABLED    = 0;
  $SiteDefs::ENSEMBL_AC_ENABLED     = 0;

  #$SiteDefs::ENSEMBL_TOOLS_JOB_DISPATCHER = { 'Blast' => 'NcbiBlast' }; #in public-plugins/tools/conf, defined 'Blast' => ''
  #$SiteDefs::NCBIBLAST_REST_ENDPOINT  = 'http://www.ebi.ac.uk/Tools/services/rest/ncbiblast'; 


  push @SiteDefs::ENSEMBL_HTDOCS_DIRS, # Needed due to EG plugin
    $SiteDefs::ENSEMBL_SERVERROOT.'/biomart-perl/htdocs';


  #----------
  #GeoLiteCity database file
  $SiteDefs::GEOCITY_DAT = $SiteDefs::ENSEMBL_SERVERROOT.'/geocity/GeoLiteCity.dat';


  #----------
  # Species stuff

    $SiteDefs::ENSEMBL_DATASETS = [sort qw(
      Arabidopsis_thaliana
      Brachypodium_distachyon
      Oryza_brachyantha
      Oryza_glaberrima
      Oryza_indica
      Oryza_indicair8
      Oryza_sativa
      Oryza_aus
      Oryza_carolina
      Oryza_barthii
      Oryza_glumaepatula
      Oryza_meridionalis
      Oryza_nivara
      Oryza_punctata
      Oryza_longistaminata
      Sorghum_bicolor
      Oryza_rufipogon
      Leersia_perrieri
    ), qw(
    )];

 $SiteDefs::PRODUCTION_NAMES = [sort qw(
        arabidopsis_thaliana
        brachypodium_distachyon
        leersia_perrieri
        oryza_barthii
        oryza_brachyantha
        oryza_glaberrima
        oryza_glumaepatula
        oryza_indica
	oryza_indicair8
        oryza_longistaminata
        oryza_meridionalis
        oryza_nivara
        oryza_punctata
        oryza_rufipogon
        oryza_sativa
	oryza_carolina
	oryza_aus
        sorghum_bicolor
    )];

%SiteDefs::__species_aliases         = (
	'arabidopsis_thaliana' => [ 'arabidopsis_thaliana' ],
        'brachypodium_distachyon' => [ 'brachypodium_distachyon' ],
        'leersia_perrieri' => [ 'leersia_perrieri' ],
        'oryza_barthii' => [ 'oryza_barthii' ],
        'oryza_brachyantha' => [ 'oryza_brachyantha' ],
        'oryza_glaberrima' => [ 'oryza_glaberrima' ],
        'oryza_glumaepatula' => [ 'oryza_glumaepatula' ],
        'oryza_indica' => ['oryza_indica'],
        'oryza_indicair8' => ['oryza_indicair8'],
        'oryza_longistaminata' => ['oryza_longistaminata'],
        'oryza_meridionalis' => ['oryza_meridionalis'],
        'oryza_nivara' => ['oryza_nivara'],
        'oryza_punctata' => ['oryza_punctata'],
        'oryza_rufipogon' => ['oryza_rufipogon'],
        'oryza_sativa' => ['oryza_sativa'],
        'oryza_carolina' => ['oryza_carolina'],
        'oryza_aus' => ['oryza_aus'],
        'sorghum_bicolor' => ['sorghum_bicolor'],
    );

  $SiteDefs::ENSEMBL_PRIMARY_SPECIES    = 'oryza_sativa'; # Default
  $SiteDefs::ENSEMBL_SECONDARY_SPECIES  = 'arabidopsis_thaliana';

}

1;
