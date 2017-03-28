#!/usr/local/bin/perl

# script to dump out all genes, transcripts, exons, introns, cdses from specified ensembl databases.
# specify the db connections with options, default to cabot
# Specify the output directory with --output_dir=/path/to/output. Defaults to /tmp.
# input database name from command line arguments

=pod

=head1 NAME

sharon_gff_dump.pl - dump out all genes, transcripts, exons, introns, cdses from specified ensembl databases.

=head1 SYNOPSIS

   sharon_gff_dump.pl [options] dbname1 dbname2 ...

Options:

  -h|--help             Show brief help and exit.
  -d|--dbhost		database host, default to cabot
  -u|--user             database useri, default to gramene_web
  --pass          	database password, default to the one for gramene_web
  --port		database port, default to 3306       
  -o|--output_dir       path to output directory, default to /tmp

=cut


use lib map { "/usr/local/ensembl-live/$_" } 
	qw ( 	bioperl-live 
		modules 
		ensembl/modules 
		ensembl-external/modules 
		ensembl-draw/modules 
		ensembl-compara/modules 
		ensembl-variation/modules/ 
		gramene-live/ensembl-plugins/gramene/modules 
		gramene-live/ensembl-plugins/maize/modules 
		conf);

use FindBin qw($Bin);
use lib "$Bin";

use strict;
use warnings;
use ExportView::GFF3Exporter;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::Registry;

my %args = (
	'output_dir' => '/tmp',
	'user' => 'gramene_web',
	'pass' => 'gram3n3',
	'port' => 3306,
	'dbhost' => 'cabot',	
);

GetOptions(
	\%args,
	'output_dir=s',
	'user=s',
	'pass=s',
	'dbhost=s',
	'port=i',	
	'help',
);

my @db_adaptors = @ARGV;

pod2usage( -verbose => 2) if $args{'help'};

print STDERR "Sending files to $args{'output_dir'}\n\n";

foreach my $dbname (@db_adaptors) {

    my $key = $dbname;

	my $tdba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
                               -user   => $args{user},
                               -pass   => $args{pass},
                               -dbname => $dbname,
                               -host   => $args{dbhost},
                               -port   => $args{port},
                               -driver => 'mysql',
                               );

    my $adaptor = $tdba->get_SliceAdaptor;

	my $slices = $adaptor->fetch_all('toplevel');

	my $output_file = $args{'output_dir'} . '/' . $key . ".gff";

	print STDERR "\tDumping $key ($dbname) to $output_file (" , scalar(@$slices), " slices)\n";

	open (my $gff, '>' . $output_file);

    eval {
	    my $exporter = ExportView::GFF3Exporter->new('debug' => 1);
	    $exporter->header($gff, $adaptor->db->get_MetaContainer->get_common_name(), $adaptor->db->get_MetaContainer->get_genebuild());
    	$exporter->export_genes_from_slices($gff, @$slices);
    };

    close $gff;

    if ($@) {
    	unlink $output_file;
    	print STDERR "BLAMMO! DUMP FAILED: $@. Deleting log!\n";
	exit;
    }

}

print "Success!\n";
