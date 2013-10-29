#!/usr/local/bin/perl

use strict;
use warnings;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

my $dbname = $ARGV[0]; # name of database
my $file = $ARGV[1];   # file containing genes to replace stops for
my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -user => 'ensembl_rw',
    -host => 'cabot',
    -port => '3306',    
    -pass => '()ryz@',
    -dbname => $dbname,
    -db_version=>64);

# open file
open(my $fh, '<', $file) or die("Unable to open file $file");
print $dba->dbc()->dbname();
# read file into an array
my @gene_ids = <$fh>;

# close file 
close($fh);
my $aa = $dba->get_AttributeAdaptor();
for my $gene_id (@gene_ids) {
    chomp ($gene_id);
    my $gene = $dba->get_GeneAdaptor()->fetch_by_stable_id($gene_id);
    if($gene) {
    print "Checking gene ".$gene->dbID."/".$gene->stable_id()."\n";
    for my $transcript (@{ $gene->get_all_Transcripts()}) {

        my $translation = $transcript->translation();
        my $idx = index($translation->seq(),'*',0);
        while($idx!=-1) {
                my $pidx = $idx+1;
		$aa->store_on_Translation($translation,
			      [
			       Bio::EnsEMBL::Attribute->new(
				   -CODE        => "amino_acid_sub",
				   -NAME        => "Amino acid substitution",
				   -DESCRIPTION => "Some translations have been manually curated for amino acid substitiutions. For example a stop codon may be changed to an amino acid in order to prevent premature truncation, or one amino acid can be substituted for another.",
				   -VALUE => "$pidx $pidx X")]);
         	$idx = index($translation->seq(),'*',$idx+1);
	}
    }
    } else {
        print "Could not find gene $gene_id\n";
    }
}

