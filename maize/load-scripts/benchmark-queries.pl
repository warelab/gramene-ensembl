#!/usr/local/bin/perl

use strict;
use warnings;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::MiscFeature;
use Bio::EnsEMBL::Attribute;
use Bio::EnsEMBL::SimpleFeature;
use Bio::EnsEMBL::MiscSet;

use FindBin qw($Bin);

use Benchmark qw(:all);

use vars qw($ENS_DBA);

MAIN: {
    my $registry_file = "$Bin/benchmark.registry";
    if (! -e $registry_file) {
        die "File does not exist: $registry_file\n";
    }
    Bio::EnsEMBL::Registry->load_all($registry_file);

    routine($ARGV[0]);
}

sub routine {
    my ($species) = @_;
    $ENS_DBA = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'core');
    # print "\nDATABASE VERSION $species\n";
    # print '=' x 60, "\n";    
    my $misc_feature_adaptor = $ENS_DBA->get_adaptor('MiscFeature');
    my $attribute_adaptor    = $ENS_DBA->get_adaptor('Attribute');

    my @features = @{ $misc_feature_adaptor->fetch_all_by_attribute_type_value('name', 'c0003D12') };
    # print "Got ", scalar @features, " features\n";
}

1;

