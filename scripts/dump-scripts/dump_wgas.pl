#!/usr/local/bin/perl -w


BEGIN {
    $ENV{'GrameneDir'} ||= '/usr/local/gramene/';
    $ENV{'HOME'} ||= '/usr/local/';
}


use lib map { $ENV{HOME} . "/$_" }
            qw ( 
		 ensembl-live/ensembl/modules
		 ensembl-live/ensembl-compara/modules
		 bioperl-live
		 gramene-ensembl/gramene/modules
		);


#GFF3ExporterWeix.pm in gramene-svn-2/ensembl-plugins/maize/modules
        


use strict;
use Bio::EnsEMBL::Registry;
use Bio::AlignIO;
use Bio::EnsEMBL::Utils::Exception qw(throw);

use Getopt::Long;
use Pod::Usage;

use ExportView::GFF3ExporterWeix;

=head1 NAME
 
   gene_overlap_wga.pl

=head1 SYNOPSIS

    perl gene_overlap_wga.pl [OPTIONS]


=head1 DESCRIPTION

    To answer the question:
    1.How many genes in each genome are touched by an HSP in the chainNet alignments and overlap in coding region for x%(--alignment_type BLASTZ_NET)


    perl gene_overlap_HSP.pl --reg_conf /usr/local/gramene_ensembl/conf/ensembl.registry --query_species Arabidopsis_thaliana --target_species Oryza_sativa_japonica --set_of_species "Oryza_sativa_japonica:Arabidopsis_thaliana" --dbname compara --alignment_type BLASTZ_NET

    perl gene_overlap_HSP.pl --reg_conf ~apurva/gramene_comparaBH/conf/ensembl.registry --query_species Zea_mays_accel -target_species Oryza_sativa_japonica --set_of_species "Zea_mays_accel:Oryza_sativa_japonica" --dbname os_accel_compara --alignment_type SYNTENY --q_coord SUPERCONTIG --t_gene_src tigr

=head1 OPTIONS

    --help           print help information
    --man            print the documentaion
    --reg_conf       ensembl registry file
    --dbname         the designated species name for compara database in the registry file
    --query_species  the designated species name for query database in the registry file
    --target_species the designated species name for target database in the registry file
    --set_of_species the set of species seperated by ":"

    --alignment_type SYNTENY or BLASTZ_NET
    --q_coord        the query species coordinate system, default is "chromosome"
    --t_coord        the target species coordinate system, default is "chromosome"
    --q_gene_src     the query species gene set source, default is all the genes in the database
    --t_gene_src     the target species gene set source, default is all the genes in the database
    --q_gene_biotype the query species gene biotype, if specified only genes of this biotype will be included in the analysis 
    --t_gene_biotype the target species gene biotype, if specified only genes of this biotype will be included in the analysis 

    --perc          the overlaping coding region percentage

=head1 AUTHORS

    Sharon Wei

=cut

my ($help, $man);
my $reg_conf;
my $set_of_species;
my $dbname;
#my ($query_species, $target_species, $q_coord, $t_coord, 
#    $q_gene_src, $t_gene_src, $q_gene_biotype, $t_gene_biotype);
my $species;
my $target_species;
my $coord;
my $gene_src;
my $gene_biotype;
my $alignment_type;
my $perc;
my $DEBUG;


GetOptions(
	   "help" => \$help,
	   "man" => \$man,
	   "reg_conf=s" => \$reg_conf,
	   "dbname=s" => \$dbname, 
	   "species=s" => \$species,
	   "set_of_species=s" => \$set_of_species,
	   "alignment_type=s" => \$alignment_type,
	   "coord=s" => \$coord,
	   "debug"           => \$DEBUG,

  );

pod2usage( -exitval=>2, -verbose=>2 ) if ($help || $man);

my $reg = "Bio::EnsEMBL::Registry";

#print "debug=$DEBUG\n";
$reg->no_version_check(1);
$reg->load_all($reg_conf);

my $compara_DBAdaptor = $reg->get_DBAdaptor($dbname, 'compara');

my $database_name = $compara_DBAdaptor->dbc->dbname;

my $genome_db_adaptor = $compara_DBAdaptor->get_GenomeDBAdaptor();
throw("cannot create adaptor for genome db connecting to <$dbname>")
    if (!$genome_db_adaptor);

my $dnafrag_adaptor = $compara_DBAdaptor->get_DnafragAdaptor();
throw("cannot create adaptor for dnafrag connecting to <$dbname>")
    if (!$dnafrag_adaptor);

my $mlss_adaptor = $compara_DBAdaptor->get_MethodLinkSpeciesSetAdaptor();
throw("cannot create adaptor for mlss connecting to <$dbname>")
    if (!$mlss_adaptor);

my $gab_adaptor = $compara_DBAdaptor->get_GenomicAlignBlockAdaptor();
throw("cannot create adaptor for GAB for connecting to <$dbname>")
    if (!$gab_adaptor);

my $DBAdaptor = $reg->get_DBAdaptor($species, "core");

my $slice_adaptor = $DBAdaptor->get_SliceAdaptor;

my $genome_dbs;
my $reference_species = lc $species;
$reference_species =~ s/\s+/_/g;
my $genomes_info;
$genomes_info->{dbname} = $database_name;

foreach my $species (split(":", $set_of_species)) {
    
    my $uniform_species = lc $species;
    $uniform_species =~ s/\s+/_/g;
    
   if($uniform_species eq $reference_species){
       $genomes_info->{ref} = $uniform_species;
   }else{
       $genomes_info->{nonref} = $uniform_species;
   }

    my $genome_db = $genome_db_adaptor->fetch_by_registry_name($species);
    
    #print "genomeb_db_id for $species is ", $genome_db->dbID;
   # Add Bio::EnsEMBL::Compara::GenomeDB object to the list
    push(@$genome_dbs, $genome_db);

    my $seq_listref = $dnafrag_adaptor->fetch_all_by_GenomeDB_region($genome_db);
    $genomes_info->{genomes}->{$uniform_species}->{dnafrag} = $seq_listref;

}
#my $num_species = scalar(@$genome_dbs);

warn ("$alignment_type, ". (join ",", map {$_->name} @{$genome_dbs}). "\n");
my $mlss = $mlss_adaptor->fetch_by_method_link_type_GenomeDBs($alignment_type, $genome_dbs);

print "mlssid is ", $mlss->dbID ; 
  
#print"chr\tTotalGenes\tHSPGenes\tHSPGenesCoding$perc\n";
#get_genes_overlapHSP($query_species, $q_coord, $q_gene_src, $q_gene_biotype);

#print_header_info();

dump_wga_gff3($slice_adaptor, $coord, $gab_adaptor, $mlss, $genomes_info);



sub dump_wga_gff3{

  my $slice_adaptor = shift;
  my $coord = shift || 'toplevel';
  my $gab_adaptor = shift;
  my $mlss = shift;
  my $genomes_info = shift;

  #print $mlss->dbID;
  my $gff3Exporter = ExportView::GFF3ExporterWeix->new();

  my @slices = @{$slice_adaptor->fetch_all($coord)};

  $gff3Exporter->wga_header(\*STDOUT, $species, $alignment_type, $genomes_info);

  foreach my $slice (@slices) {

      
      $gff3Exporter->export_pairwise_wga_from_slice(\*STDOUT, $slice, $gab_adaptor, $mlss);
 
  }
   
}



