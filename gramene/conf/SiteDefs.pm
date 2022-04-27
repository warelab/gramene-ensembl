
package EnsEMBL::Gramene::SiteDefs;
use strict;

sub update_conf {
    $SiteDefs::ENSEMBL_PORT           = 8888;
    $SiteDefs::ENSEMBL_SERVERNAME     = 'oryza-ensembl.gramene.org';

    $SiteDefs::ENSEMBL_PRIMARY_SPECIES   = 'Oryza_sativa';
    $SiteDefs::ENSEMBL_SECONDARY_SPECIES = 'Arabidopsis_thaliana';

    $SiteDefs::ENSEMBL_BASE_URL          =  $SiteDefs::ENSEMBL_SERVERNAME;
    $SiteDefs::SITE_RELEASE_VERSION = 3.1;
    $SiteDefs::SITE_RELEASE_DATE    = 'Mar 2022';
	#'Nov 2021';


    $SiteDefs::EG_DIVISION      = 'plants';
    $SiteDefs::ENSEMBL_SITETYPE = 'GrameneOryza';
    $SiteDefs::SITE_NAME        = 'GrameneOryza';
    $SiteDefs::SITE_FTP         = 'http://ftp.gramene.org/oge';

    $SiteDefs::ENSEMBL_USER       = 'nobody';#getpwuid($>);
    $SiteDefs::ENSEMBL_GROUP      = 'nobody';#getgrgid($));
    $SiteDefs::ENSEMBL_SERVERADMIN            = 'weix@cshl.edu';
    $SiteDefs::ENSEMBL_MAIL_SERVER       = 'localhost';

    @SiteDefs::ENSEMBL_PERL_DIRS =
      ( $SiteDefs::ENSEMBL_WEBROOT.'/perl',
        $SiteDefs::ENSEMBL_SERVERROOT.'/eg-web-common/perl',
        $SiteDefs::ENSEMBL_SERVERROOT.'/eg-web-plants/perl',
      );


    $SiteDefs::ENSEMBL_DATASETS = [sort qw(
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
      Oryza_rufipogon
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

     $SiteDefs::NCBIBLAST_REST_ENDPOINT = 'http://brie:5202';

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

 our $VCF_LIB = $SiteDefs::ENSEMBL_SERVERROOT.'/vcftools/perl';
 push @SiteDefs::ENSEMBL_LIB_DIRS, $VCF_LIB;

 $SiteDefs::EBEYE_REST_ENDPOINT = 'https://data.gramene.org/oryza-ebeye3';

}

1;
