package EnsEMBL::Gramene::SiteDefs;
use strict;

# These are the gramene-specific edits to the main Ensembl SiteDefs.pm file
sub update_conf {

  $SiteDefs::ENSEMBL_SERVERNAME             = 'maize-pangenome-ensembl.gramene.org';
#  $SiteDefs::ENSEMBL_SERVERNAME             = 'dev.gramene.org';

  $SiteDefs::ENSEMBL_BASE_URL     = $SiteDefs::ENSEMBL_SERVERNAME;
  $SiteDefs::SITE_RELEASE_VERSION = 2; 
  $SiteDefs::SITE_RELEASE_VERSION_EG = 32;
  #$SiteDefs::ENSEMBL_VERSION = 65;
  $SiteDefs::SITE_RELEASE_DATE    = 'Jan 2019';
  $SiteDefs::SITE_NAME            = 'Gramene';
  $SiteDefs::SITE_FTP             = 'ftp://ftp.gramene.org/pub';
  $SiteDefs::GRAMENE_FTP_URL	  = 'ftp://ftp.gramene.org/pub';
#  $SiteDefs::OGE_FTP_URL          = 'ftp://ftp.gramene.org/pub/gramene/oge/release-current'
  $SiteDefs::PE_URL             = 'http://plants.ensembl.org';
  $SiteDefs::ENSEMBL_PORT       = 8888;
  $SiteDefs::ENSEMBL_PROXY_PORT = 80; # Port used for self-referential URLs
  $SiteDefs::ENSEMBL_USER       = 'nobody';#getpwuid($>);          
  $SiteDefs::ENSEMBL_GROUP      = 'nobody';#getgrgid($));           

  $SiteDefs::ENSEMBL_SERVERADMIN            = 'weix@cshl.edu';
  $SiteDefs::ENSEMBL_MAIL_SERVER            = 'localhost';

  $SiteDefs::SAMTOOLS_DIR = $SiteDefs::ENSEMBL_SERVERROOT.'/samtools'; 

  $SiteDefs::ENSEMBL_DEBUG_FLAGS             = 0; # 24;
  $SiteDefs::ENSEMBL_LONGPROCESS_MINTIME     = 10;

#  $SiteDefs::EBEYE_REST_ENDPOINT     = 'http://brie:12051';
  $SiteDefs::EBEYE_REST_ENDPOINT = 'http://data.gramene.org/panmaize-ebeye';

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

  #----------
  # User database
  $SiteDefs::ENSEMBL_USERDB_NAME = 'ensembl_accounts';  #changed to ensembl_accounts lately
  #'ensembl_web_user_db_31_57';
  $SiteDefs::ENSEMBL_USERDB_USER = 'gramene_web';
  $SiteDefs::ENSEMBL_USERDB_HOST = 'cabot.cshl.edu';
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
  $SiteDefs::ENSEMBL_BLAST_ENABLED = 0; # Creates header link for blast
  $SiteDefs::ENSEMBL_MART_ENABLED = 0; # And mart

  $SiteDefs::ENSEMBL_VEP_ENABLED    = 0;
  $SiteDefs::ENSEMBL_AC_ENABLED     = 0;

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
  push @{$SiteDefs::ENSEMBL_DATASETS}, 'Zea_maysw22';
  push @{$SiteDefs::ENSEMBL_DATASETS}, 'Zea_maysph207';
  push @{$SiteDefs::ENSEMBL_DATASETS}, 'Zea_maysb73';
  push @{$SiteDefs::ENSEMBL_DATASETS}, 'Zea_mayscml247';
  push @{$SiteDefs::ENSEMBL_DATASETS}, 'Zea_mayscml333';
  push @{$SiteDefs::ENSEMBL_DATASETS}, 'Zea_mayshp301';
  push @{$SiteDefs::ENSEMBL_DATASETS}, 'Zea_maysil14h';
  push @{$SiteDefs::ENSEMBL_DATASETS}, 'Zea_maysky21';
  push @{$SiteDefs::ENSEMBL_DATASETS}, 'Setaria_italica';
  push @{$SiteDefs::ENSEMBL_DATASETS}, 'Glycine_max';
  push @{$SiteDefs::ENSEMBL_DATASETS}, 'Brachypodium_distachyon';
  push @{$SiteDefs::ENSEMBL_DATASETS}, 'Sorghum_bicolor';
  push @{$SiteDefs::ENSEMBL_DATASETS}, 'Vitis_vinifera';
  push @{$SiteDefs::ENSEMBL_DATASETS}, 'Physcomitrella_patens';
  push @{$SiteDefs::ENSEMBL_DATASETS}, 'Chlamydomonas_reinhardtii';
  push @{$SiteDefs::ENSEMBL_DATASETS}, 'Selaginella_moellendorffii';

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
