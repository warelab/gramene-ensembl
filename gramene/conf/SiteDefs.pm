package EnsEMBL::Gramene::SiteDefs;

use strict;

# These are the gramene-specific edits to the main Ensembl SiteDefs.pm file

sub update_conf {

  #$SiteDefs::ENSEMBL_SERVERNAME             = 'ensembl.sorghumbase.org';
  $SiteDefs::ENSEMBL_SERVERNAME             = 'ensembl-dev.sorghumbase.org';
  
  $SiteDefs::GRM_SERVERNAME = 'sorghumbase.org';

  $SiteDefs::ENSEMBL_MAX_PROCESS_SIZE     = 4000000; 
  $SiteDefs::ENSEMBL_BASE_URL     = $SiteDefs::ENSEMBL_SERVERNAME;
#  $SiteDefs::ENSEMBL_STATIC_SERVER  = $SiteDefs::ENSEMBL_SERVERNAME; 
	#the ensembl-webcode/modules/Image/Minifier.pm in v87 uses ENSEMBL_STATIC_SERVER, but it was not defined anywhere, maybe a bug in that module, define it here
  $SiteDefs::SITE_RELEASE_VERSION = 7;  #this is sorghumbase version1, but the outgroup db cores are _2_87
  $SiteDefs::SITE_RELEASE_VERSION_EG = 43;
  $SiteDefs::SITE_RELEASE_DATE    = 'Sep 2023';

  $SiteDefs::LARGE_SPECIES_SET = 1;
  $SiteDefs::SITE_NAME            = 'Sorghumbase';
  $SiteDefs::ENSEMBL_SITETYPE = 'Sorghumbase';
	#$SiteDefs::ENSEMBL_SITETYPE = 'Ensembl Plants';
  $SiteDefs::SITE_FTP             = 'https://ftp.sorghumbase.org';
  $SiteDefs::PE_URL             = 'https://'.$SiteDefs::ENSEMBL_SERVERNAME;
  $SiteDefs::ENSEMBL_PORT       = 88;
  $SiteDefs::ENSEMBL_PROXY_PORT = 80; # Port used for self-referential URLs
  $SiteDefs::ENSEMBL_USER       = 'nobody';#getpwuid($>);          
  $SiteDefs::ENSEMBL_GROUP      = 'nobody';#getgrgid($));           

  $SiteDefs::ENSEMBL_SERVERADMIN            = 'weix@cshl.edu';
  $SiteDefs::ENSEMBL_MAIL_SERVER            = 'localhost';

  $SiteDefs::MAX_PROCESS_SIZE = 2000000;

  $SiteDefs::VCFTOOLS_PERL_LIB    = $SiteDefs::ENSEMBL_SERVERROOT . "/vcftools/src/perl/";

  $SiteDefs::GRAPHIC_TTF_PATH = '/usr/share/fonts/msttcore/';

  $SiteDefs::SPECIES_IMAGE_DIR  = $SiteDefs::ENSEMBL_SERVERROOT . '/gramene-live/gramene/htdocs/i/species/48/';

# start of copy of eg plants

$SiteDefs::EG_DIVISION      = 'plants';

    @SiteDefs::ENSEMBL_PERL_DIRS =
      ( $SiteDefs::ENSEMBL_WEBROOT.'/perl',
        $SiteDefs::ENSEMBL_SERVERROOT.'/eg-web-common/perl',
        $SiteDefs::ENSEMBL_SERVERROOT.'/eg-web-plants/perl',
      );


    $SiteDefs::DOCSEARCH_INDEX_DIR = $SiteDefs::ENSEMBL_SERVERROOT.'/eg-web-plants/data/docsearch';

    $SiteDefs::ENA_COLLECTION_ID = 224;

    $SiteDefs::ENA_SAMPLE_SEQ = "MDDCRFETSELQASVMISTPLFTDSWSSCNTANCNGSIKIHDIAGITYVAIPAVSMIQLGNLVGLPVTGDVLFPGLSSDEPLPMVDAAILKLFLQLKIKEGLELELLGKKLVVITGHSTGGALAAFTALWLLSQSSPPSFRVFCITFGSPLLGNQSLSTSISRSRLAHNFCHVVSIHDLVPRSSNEQFWPFGTYLFCSDKGGVCLDNAGSVRLMFNILNTTATQNTEEHQRYGHYVFTLSHMFLKSRSFLGGSIPDNSYQAGVALAVEALGFSNDDTSGVLVKECIETATRIVRAPILRSAELANELASVLPARLEIQWYKDRCDASEEQLGYYDFFKRYSLKRDFKVNMSRIRLAKFWDTVIKMVETNELPFDFHLGKKWIYASQFYQLLAEPLDIANFYKNRDIKTGGHYLEGNRPKRYEVIDKWQKGVKVPEECVRSRYASTTQDTCFWAKLEQAKEWLDEARKESSDPQRRSLLREKIVPFESYANTLVTKKEVSLDVKAKNSSYSVWEANLKEFKCKMGYENEIEMVVDESDAMET";

    $SiteDefs::GXA = 1;

    $SiteDefs::ENSEMBL_HMMER_ENABLED  = 1;

#end of eg-plant conf

  $SiteDefs::SAMTOOLS_DIR = $SiteDefs::ENSEMBL_SERVERROOT.'/samtools'; 

  $SiteDefs::ENSEMBL_DEBUG_FLAGS             = 24;
  $SiteDefs::ENSEMBL_LONGPROCESS_MINTIME     = 10;


  $SiteDefs::NCBIBLAST_REST_ENDPOINT = 'http://brie:5202';
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
    'host'        => 'cabot',                                                                       # Database host (defaults to ensembldb.ensembl.org)
    'user'        => 'weix',                                                                       # Defaults to 'anonymous'
    'password'    => 'warelab',                                                                       # Not used by default
    'port'        => 3306,                                                                       # Defaults to 5306
    'fork'        => 4,                                                                           # Enable forking, using 4 forks
  };

  #----------
  # User database
  $SiteDefs::ENSEMBL_USERDB_NAME = 'ensembl_accounts';  #changed to ensembl_accounts lately
  $SiteDefs::ENSEMBL_USERDB_USER = 'weix';
  $SiteDefs::ENSEMBL_USERDB_HOST = 'cabot.cshl.edu';
  $SiteDefs::ENSEMBL_USERDB_PORT =  3306;
  $SiteDefs::ENSEMBL_USERDB_PASS = 'warelab';

  #----------
  # Logging
  $SiteDefs::ENSEMBL_LOGDIR = $SiteDefs::ENSEMBL_SERVERROOT."/logs";
  $SiteDefs::ENSEMBL_PIDFILE   = "$SiteDefs::ENSEMBL_LOGDIR/httpd.pid";
  $SiteDefs::ENSEMBL_ERRORLOG  = "$SiteDefs::ENSEMBL_LOGDIR/error.log";
  $SiteDefs::ENSEMBL_CUSTOMLOG = "$SiteDefs::ENSEMBL_LOGDIR/access.log combined";

  #----------
  # Mart/Blast
  $SiteDefs::ENSEMBL_BLAST_ENABLED = 1; # Creates header link for blast
  $SiteDefs::ENSEMBL_BLAST_BY_SEQID = 1; # blast on gene page sequence
  $SiteDefs::ENSEMBL_MART_ENABLED = 0; # And mart

  $SiteDefs::ENSEMBL_VEP_ENABLED    = 0;
  $SiteDefs::ENSEMBL_AC_ENABLED     = 0;
  $SiteDefs::ENSEMBL_IDM_ENABLED    = 0;
  $SiteDefs::ENSEMBL_FC_ENABLED     = 0;
  $SiteDefs::ENSEMBL_LD_ENABLED     = 0;
  $SiteDefs::ENSEMBL_VP_ENABLED     = 0;
  $SiteDefs::ENSEMBL_DS_ENABLED     = 0;

  $SiteDefs::NCBIBLAST_REST_ENDPOINT = 'http://squam:2502'; #brie:5202';
#'http://brie:5202';
#'http://www.ebi.ac.uk/Tools/services/rest/ncbiblast';
  $SiteDefs::EBEYE_REST_ENDPOINT     =  'http://localhost:9806'; #'http://brie8:9804';
#'https://data.sorghumbase.org/ebeye' . $SiteDefs::SITE_RELEASE_VERSION;
#'http://www.ebi.ac.uk/ebisearch/ws/rest';
#$SiteDefs::EBEYE_REST_ENDPOINT     = 'http://data.gramene.org/sorghum-ebeye' . $SiteDefs::SITE_RELEASE_VERSION;


  #----------
  #GeoLiteCity database file
  $SiteDefs::GEOCITY_DAT = $SiteDefs::ENSEMBL_SERVERROOT.'/geocity/GeoLiteCity.dat';


  #----------
  # Species stuff

  $SiteDefs::ENSEMBL_PRIMARY_SPECIES    = 'Oryza_sativa'; # Default
  $SiteDefs::ENSEMBL_SECONDARY_SPECIES  = 'Arabidopsis_thaliana';

$SiteDefs::PRODUCTION_NAMES = [sort qw(
      arabidopsis_thaliana
      chlamydomonas_reinhardtii
      oryza_sativa
      selaginella_moellendorffii
      sorghum_bicolor
      sorghum_rio
      sorghum_tx2783pac
      sorghum_tx430nano
      sorghum_tx436pac
      vitis_vinifera
      zea_maysb73v4
    ),
        qw(
        sorghum_austrcf317961
        sorghum_is12661
        sorghum_is19953
        sorghum_is36143
        sorghum_is8525
        sorghum_is929
        sorghum_pi525695
        sorghum_pi532566
        sorghum_pi536008
        sorghum_r93194522
        sorghum_s3691
        zea_maysb73
    ),
        qw(
        sorghum_ji2731
        sorghum_353
        populus_trichocarpa
    ),
        qw(
        sorghum_chineseamber
        sorghum_grassl
        sorghum_leoti
        sorghum_pi229841
        sorghum_pi297155
        sorghum_pi329311
        sorghum_pi506069
        sorghum_pi510757
        sorghum_pi655972
        sorghum_riouncc 
    ),

        qw(
        sorghum_bicolorv5
    )
    ];

  $SiteDefs::ENSEMBL_TOOLS_LIST = [
    'Blast'             => 'BLAST/BLAT',
    #'VEP'               => 'Variant Effect Predictor',
    #'FileChameleon'     => 'File Chameleon',
    #'AssemblyConverter' => 'Assembly Converter',
    #'IDMapper'          => 'ID History Converter',
    #'AlleleFrequency'   => 'Allele Frequency Calculator',
    #'VcftoPed'          => 'VCF to PED Converter',
    #'DataSlicer'        => 'Data Slicer',
    #'VariationPattern'  => 'Variation Pattern Finder',
    #'LD'                => 'Linkage Disequilibrium Calculator',
  ];


# The actually configuration for the expressionAtlas restful API end point 'GXA_REST_URL' is hard-coded
# in public-plugins/widgets/htdocs/widgets/95_GXA.js 
# Not in the public-plugins/widgets/conf/SiteDefs.pm 
#  $SiteDefs::GXA_REST_URL = 'https://www.ebi.ac.uk/gxa/json/expressionData?geneId=';#'http://wwwdev.ebi.ac.uk/gxa/json/expressionData?geneId=';
#  $SiteDefs::GXA_EBI_URL  = 'https://www.ebi.ac.uk/gxa/resources';#'http://wwwdev.ebi.ac.uk/gxa/resources'; #dev  environment for GXA for pre testing their release


}

1;
