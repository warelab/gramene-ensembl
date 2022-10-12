
package EnsEMBL::Gramene::SiteDefs;
use strict;

sub update_conf {
    $SiteDefs::ENSEMBL_PORT           = 8888;
    #$SiteDefs::ENSEMBL_PROXY_PORT     = 80;
    $SiteDefs::ENSEMBL_SERVERNAME     = 'oryza-ensembl.gramene.org';
    #$SiteDefs::ENSEMBL_PROTOCOL          = 'https';
    $SiteDefs::ENSEMBL_BASE_URL          =  'https://oryza-ensembl.gramene.org/';
    $SiteDefs::SITE_RELEASE_VERSION = 5;
    $SiteDefs::SITE_RELEASE_DATE    = 'Oct 2022';
	#'May 2022';
	#'Nov 2021';


$SiteDefs::DIVISION         = 'plants';
    $SiteDefs::EG_DIVISION      = 'plants';
$SiteDefs::SUBDOMAIN_DIR    = 'plants';
    $SiteDefs::ENSEMBL_SITETYPE = 'Ensembl Plants';
    $SiteDefs::SITE_NAME        = 'GrameneOryza';
    $SiteDefs::SITE_FTP         = 'http://ftp.gramene.org/oge';
 $SiteDefs::PE_URL             = 'http://plants.ensembl.org';

    $SiteDefs::ENSEMBL_USER       = 'nobody';#getpwuid($>);
    $SiteDefs::ENSEMBL_GROUP      = 'nobody';#getgrgid($));
    $SiteDefs::ENSEMBL_SERVERADMIN            = 'weix@cshl.edu';
    $SiteDefs::ENSEMBL_MAIL_SERVER       = 'localhost';

 $SiteDefs::MAX_PROCESS_SIZE = 2000000;
 # this is defined in default.ini
 #$SiteDefs::GRAPHIC_TTF_PATH     = "/usr/share/fonts/msttcorefonts/";


  $SiteDefs::SAMTOOLS_DIR = $SiteDefs::ENSEMBL_SERVERROOT.'/samtools';

  $SiteDefs::ENSEMBL_DEBUG_FLAGS             = 0; # 24;
  $SiteDefs::ENSEMBL_LONGPROCESS_MINTIME     = 10;


  $SiteDefs::DATAFILE_ROOT        = '/usr/local';                                  ## Base path for ro data files
  $SiteDefs::DATAFILE_BASE_PATH   = '/usr/local/vcf';

$SiteDefs::EBEYE_REST_ENDPOINT = 'https://data.gramene.org/oryza-ebeye4';

# This should be set by the plugin, but not working. Need to hard code!
    #@SiteDefs::ENSEMBL_PERL_DIRS =
    #  ( $SiteDefs::ENSEMBL_WEBROOT.'/perl',
    #    $SiteDefs::ENSEMBL_SERVERROOT.'/eg-web-common/perl',
    #    $SiteDefs::ENSEMBL_SERVERROOT.'/eg-web-plants/perl',
    #  );

$SiteDefs::ENSEMBL_PRIMARY_SPECIES   = 'Oryza_sativa';
$SiteDefs::ENSEMBL_SECONDARY_SPECIES = 'Arabidopsis_thaliana';

#    $SiteDefs::ENSEMBL_DATASETS 

#$SiteDefs::PRODUCTION_NAMES = ['Oryza_sativa'];

$SiteDefs::ENSEMBL_DATASETS = [sort qw(
     Oryza_aus
     Oryza_barthii
     Oryza_brachyantha
     Oryza_carolina
     Oryza_glaberrima
     Oryza_glumaepatula
     Oryza_indica
     Oryza_indicair8
     Oryza_meridionalis
     Oryza_nivara
     Oryza_punctata
     Oryza_rufipogon
     Oryza_sativa
     Oryza_sativa117425
     Oryza_sativa125619
     Oryza_sativa125827
     Oryza_sativa127518
     Oryza_sativa127564
     Oryza_sativa127652
     Oryza_sativa127742
     Oryza_sativa128077
     Oryza_sativa132278
     Oryza_sativa132424
     Oryza_sativaazucena
     Oryza_sativair64
     Oryza_sativakitaake
     Oryza_sativamh63
     Oryza_sativazs97
      Leersia_perrieri
      Arabidopsis_thaliana
      Chlamydomonas_reinhardtii
      Drosophila_melanogaster
      Selaginella_moellendorffii
      Sorghum_bicolor
      Vitis_vinifera
      Zea_maysb73
      Zea_maysb73v4
    ), qw(
    )];

    #push @SiteDefs::ENSEMBL_HTDOCS_DIRS, $SiteDefs::ENSEMBL_SERVERROOT. '/../biomarts/plants/biomart-perl/htdocs';
      
    #$SiteDefs::DOCSEARCH_INDEX_DIR = $SiteDefs::ENSEMBL_SERVERROOT.'/eg-web-plants/data/docsearch';

    #$SiteDefs::ENA_COLLECTION_ID = 224;

    #$SiteDefs::ENA_SAMPLE_SEQ = "MDDCRFETSELQASVMISTPLFTDSWSSCNTANCNGSIKIHDIAGITYVAIPAVSMIQLGNLVGLPVTGDVLFPGLSSDEPLPMVDAAILKLFLQLKIKEGLELELLGKKLVVITGHSTGGALAAFTALWLLSQSSPPSFRVFCITFGSPLLGNQSLSTSISRSRLAHNFCHVVSIHDLVPRSSNEQFWPFGTYLFCSDKGGVCLDNAGSVRLMFNILNTTATQNTEEHQRYGHYVFTLSHMFLKSRSFLGGSIPDNSYQAGVALAVEALGFSNDDTSGVLVKECIETATRIVRAPILRSAELANELASVLPARLEIQWYKDRCDASEEQLGYYDFFKRYSLKRDFKVNMSRIRLAKFWDTVIKMVETNELPFDFHLGKKWIYASQFYQLLAEPLDIANFYKNRDIKTGGHYLEGNRPKRYEVIDKWQKGVKVPEECVRSRYASTTQDTCFWAKLEQAKEWLDEARKESSDPQRRSLLREKIVPFESYANTLVTKKEVSLDVKAKNSSYSVWEANLKEFKCKMGYENEIEMVVDESDAMET";

    #$SiteDefs::GXA = 1;

    #$SiteDefs::ENSEMBL_HMMER_ENABLED  = 1;

     $SiteDefs::ENSEMBL_BLAST_ENABLED = 1; # Creates header link for blast
     $SiteDefs::ENSEMBL_BLAST_BY_SEQID = 1; # blast on gene page sequence

     $SiteDefs::NCBIBLAST_REST_ENDPOINT = 'http://squam:2502'; #brie:5202';

     $SiteDefs::ENSEMBL_TOOLS_LIST = [
    	'Blast'             => 'BLAST/BLAT',
    #	'VEP'               => 'Variant Effect Predictor',
    #'FileChameleon'     => 'File Chameleon',
    #	'AssemblyConverter' => 'Assembly Converter',
    #'IDMapper'          => 'ID History Converter',
    #'AlleleFrequency'   => 'Allele Frequency Calculator',
    #'VcftoPed'          => 'VCF to PED Converter',
    #'DataSlicer'        => 'Data Slicer',
    #'VariationPattern'  => 'Variation Pattern Finder',
    #'LD'                => 'Linkage Disequilibrium Calculator',
  	];

# our $VCF_LIB = $SiteDefs::ENSEMBL_SERVERROOT.'/vcftools/perl';
# push @SiteDefs::ENSEMBL_LIB_DIRS, $VCF_LIB;

}

1;
