#!/usr/bin/env perl

use strict;
use warnings;
use autodie;
use feature 'say';
use Getopt::Long;
use Pod::Usage;
use Readonly;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Registry;


my ( $help, $man_page);

my ($dbhost, $dbuser, $dbpass, $dbport) = qw(cabot gramene_web gram3n3 3306);
my $dbname;

my $input_gene_cnt_sql = "select count(*) from gene_member g join genome_db n using (genome_db_id)  where source_name='ENSEMBLGENE' and n.name not in ('ciona_savignyi', 'cyanidioschyzon_merolae', 'drosophila_melanogaster', 'homo_sapiens', 'saccharomyces_cerevisiae')";

#"select count(*) from gene_member where source_name='ENSEMBLGENE'";


my $tree_memeber_fam_cnt_sql = "select  count(distinct(seq_member_id)), count(distinct(root_id)), count(distinct(g.genome_db_id)) from gene_tree_node grn  join gene_tree_root gtr using (root_id) join seq_member s using (seq_member_id) join genome_db g using (genome_db_id) where tree_type='tree' and clusterset_id='default' and g.name not in ('ciona_savignyi', 'cyanidioschyzon_merolae', 'drosophila_melanogaster', 'homo_sapiens', 'saccharomyces_cerevisiae')";

#"select  count(distinct(seq_member_id)), count(distinct(root_id)) from gene_tree_node grn  join gene_tree_root gtr using (root_id) where tree_type='tree' and clusterset_id='default'";

GetOptions(
    'help' => \$help,
    'man'  => \$man_page,
    'dbhost=s' => \$dbhost,
    'dbuser=s' => \$dbuser,
    'dbpass=s' => \$dbpass,
    'dbport=s' => \$dbport,
    'dbname=s' => \$dbname,
) or pod2usage(2);

if ( $help || $man_page ) {
    pod2usage({
        -exitval => 0,
        -verbose => $man_page ? 2 : 1
    });
}; 

my $db = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(
    '-port'    => $dbport,
    '-host'    => $dbhost,
    '-user'    => $dbuser,
    '-pass'    => $dbpass,
    '-dbname'  => $dbname,
);

my $dbh = $db->dbc->db_handle;

#$dbh->do ($input_gene_cnt_sql);
my ($input_gene_cnt) = $dbh->selectrow_array ($input_gene_cnt_sql);

#$dbh->do( $tree_memeber_fam_cnt_sql );
my ($gene_tree_member_cnt, $family_cnt, $genome_cnt) =$dbh->selectrow_array($tree_memeber_fam_cnt_sql);
$dbh->disconnect;

#57,562 gene are from outgroup species: 'ciona_savignyi', 'cyanidioschyzon_merolae', 'drosophila_melanogaster', 'homo_sapiens', 'saccharomyces_cerevisiae'

printf "A total of %d GeneTree families were constructed comprising %d individual genes from %d plant genomes with %d input proteins.\n", $family_cnt, $gene_tree_member_cnt, $genome_cnt, $input_gene_cnt;




__END__

# ----------------------------------------------------

=pod

=head1 NAME

ensembl-genetree-stats.pl - cacluate gene tree statistics based on ensembl-compara database

=head1 SYNOPSIS

  ensembl-genetree-stats.pl 

Options:

  --help   Show brief help and exit
  --man    Show full documentation
  -dbuser  database connection 
  -dbpass  database connection
  -dbport  database connection
  -dbhost  database connection
  -dbname  database name

=head1 DESCRIPTION

Describe what the script does, what input it expects, what output it
creates, etc.

=head1 SEE ALSO

perl.

=head1 AUTHOR

weix E<lt>weix@cshl.eduE<gt>.

=head1 COPYRIGHT

Copyright (c) 2014 Cold Spring Harbor Laboratory

This module is free software; you can redistribute it and/or
modify it under the terms of the GPL (either version 1, or at
your option, any later version) or the Artistic License 2.0.
Refer to LICENSE for the full license text and to DISCLAIMER for
additional warranty disclaimers.

=cut
