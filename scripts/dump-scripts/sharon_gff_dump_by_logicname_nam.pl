#!/usr/local/bin/perl

# script to dump out all genes, transcripts, exons, introns, cdses from all gramene ensembl databases.
# pulls out all databases from the gramene config file that start with ensembl_ and are hosted on cabot.
# Specify the output directory with --output_dir=/path/to/output. Defaults to /tmp.

#use lib '/usr/local/gramene/lib/perl';
#use lib '/usr/local/gramene-svn-2/ensembl-plugins/maize/modules';
#use lib '/usr/local/gramene-cvs/ensembl-plugins/maize/modules/';

#use lib map { "/usr/local/ensembl-live/$_" } qw ( bioperl-live modules ensembl/modules ensembl-external/modules ensembl-draw/modules ensembl-compara/modules ensembl-variation/modules/ gramene-live/ensembl-plugins/gramene/modules gramene-live/ensembl-plugins/maize/modules conf);
use lib map { "/usr/local/ensembl-live/$_" } qw ( bioperl-live modules ensembl/modules ensembl-external/modules ensembl-draw/modules ensembl-compara/modules ensembl-variation/modules/ gramene-live/ensembl-plugins/gramene/modules gramene-live/ensembl-plugins/maize/modules conf);

use FindBin qw($Bin);
use lib "$Bin";

use strict;
use warnings;
use ExportView::GFF3ExporterNAM;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;
use Bio::EnsEMBL::Registry;

my @db_adaptors;
my ($output_dir, $logicname);

GetOptions(
	'dbname=s' => \@db_adaptors,
	'output_dir=s' => \$output_dir,
	'logicname=s' => \$logicname,
);


#qw(
#leersia_perrieri_core_38_72_03mk
#oryza_granulata_core_38_72_03mk
#oryza_longistaminata_core_38_72_03mk
#oryza_meridionalis_core_38_72_03mk
#oryza_barthii_core_38_72_03mk
#oryza_brachyantha_core_38_72_03mk
#oryza_minutabb_core_38_72_03mk
#oryza_minutacc_core_38_72_03mk
#oryza_nivara_core_38_72_03mk
#oryza_officinalis_core_38_72_03mk
#oryza_punctata_core_38_72_03mk
#oryza_rufipogon_core_38_72_03mk
#);
#leersia_perrieri_core_38_72_10mk
#oryza_barthii_core_38_72_1mk
#oryza_brachyantha_core_38_72_14mk
#oryza_glaberrima_core_38_72_1mk
#oryza_glumaepatula_core_38_72_15mk
#oryza_indica_core_38_72_2mk
#oryza_longistaminata_core_38_72_117mk
#oryza_meridionalis_core_38_72_1mk
#oryza_nivara_core_38_72_10mk
#oryza_punctata_core_38_72_12mk
#oryza_rufipogon_core_38_72_11mk
#oryza_sativa_core_38_72_7mk);

print STDERR "Sending files to $output_dir\n\n";

foreach my $dbname (@db_adaptors) {

    my $key = $dbname;

	my $tdba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
                               -user   => 'plensembl',
                               -pass   => 'AudreyII',
                               -dbname => $dbname,
                               -host   => 'bhsqldw1',
                               -port    => '3306',
                               -driver  => 'mysql',
                               );

    my $adaptor = $tdba->get_SliceAdaptor;

	my $slices = $adaptor->fetch_all('toplevel');

	my $output_file = $output_dir . '/' . $key . ".gff";

	print STDERR "\tDumping $key ($dbname) to $output_file (" , scalar(@$slices), " slices)\n";
my $gff;
	open ( $gff, '>' . $output_file) or die "Cannot open $output_file\n";

    eval {
	    my $exporter = ExportView::GFF3ExporterNAM->new('debug' => 1,
							    'source' => 'NAM');
	    $exporter->header($gff, $adaptor->db->get_MetaContainer->get_common_name(), $adaptor->db->get_MetaContainer->get_genebuild());
    	$exporter->export_genes_from_slices($gff, $slices, $logicname);
    };

    close $gff;

    if ($@) {
    	unlink $output_file;
    	print STDERR "BLAMMO! DUMP FAILED: $@. Deleting log!\n";
    }

}

print "Success!\n";
