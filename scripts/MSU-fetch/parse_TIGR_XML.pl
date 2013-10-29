#!/usr/local/bin/perl

## Contact Brian Haas (bhaas@tigr.org) with questions/comments/bugs/etc.

use lib '/usr/local/gramene/scripts/ensembl/MSU-fetch/TIGRXMLLIB/'; 
use lib '/usr/local/gramene/scripts/ensembl/MSU-fetch/haaslib/';

#use lib ("/usr/local/data/tigr/TIGRXMLLIB");
use TIGR_XML_parser;
use Gene_obj;
use strict;

my $file;
unless ($file = $ARGV[0]) {die "\n\nusage: $0 xmlfile\n\n";}

my $TIGRparser = new TIGR_XML_parser();
my $pwd = `pwd`;
chomp $pwd;
unless ($file =~ /\//) {
    #full path not specified.
    $file = "$pwd/$file";
}

$TIGRparser->capture_genes_from_assembly_xml("$file");

#print $TIGRparser->toString();

my @genes = $TIGRparser->get_genes();
my $x = 1;
foreach my $gene (@genes) {
    $gene->refine_gene_object();
    # check out the Gene_obj.pm toString() method
    # to see what simple access methods are available.
    print "$x\n" . $gene->toString();
    $x++;
}


