#!/usr/bin/env perl

use lib '/usr/local/gramene/scripts/ensembl/MSU-fetch/TIGRXMLLIB/'; 
use lib '/usr/local/gramene/scripts/ensembl/MSU-fetch/haaslib/';
#("haaslib");

use TIGR_XML_parser;
use Gene_obj;
use strict;

my $file;
unless ($file = $ARGV[0]) {die "\n\nusage: $0 xmlFile\n\n"};

my $TIGRparser = new TIGR_XML_parser();
my $pwd = `pwd`;
chomp $pwd;
$TIGRparser->capture_genes_from_assembly_xml("$pwd/$file");


my @genes = $TIGRparser->get_genes();

foreach my $gene (@genes) {
    print $gene->to_GFF3_format();
    print "\n"; #spacer between genes.
}


exit(0);
