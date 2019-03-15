#!/usr/bin/env perl



use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;

use Text::RecordParser::Tab;
use String::CamelCase qw(camelize decamelize wordsplit);

=head1 SYNOPSIS

convert_sqldump2rmt.pl  file 

 $ perl convert_sqldump2rmt.pl MySQLdumpFile 

The input file is generated from query 

SELECT g1.name as species, gm1.stable_id,  g2.name as otherspecies, gm2.stable_id as other_stable_id, hm1.perc_id, hm2.perc_id as other_perc_id, h.is_high_confidence FROM method_link_species_set INNER JOIN homology h USING(method_link_species_set_id) INNER JOIN homology_member hm1 on (h.homology_id=hm1.homology_id) INNER JOIN gene_member gm1 on (hm1.gene_member_id=gm1.gene_member_id) INNER JOIN genome_db g1 ON gm1.genome_db_id = g1.genome_db_id INNER JOIN homology_member hm2 on (h.homology_id=hm2.homology_id) INNER JOIN gene_member gm2 on (hm2.gene_member_id=gm2.gene_member_id) INNER JOIN genome_db g2 ON gm2.genome_db_id = g2.genome_db_id WHERE method_link_id = 201 AND g1.genome_db_id = 1985 AND g2.genome_db_id != 1985 order by otherspecies

and looks like

species stable_id       otherspecies    other_stable_id       perc_id other_perc_id is_high_confidence
oryza_sativa    BAC19852        aegilops_tauschii       AET0Gv20080100  27.6923 73.1707 0
oryza_sativa    BAC19852        aegilops_tauschii       AET0Gv20080400  27.3846 87.2549 0
oryza_sativa    BAC19852        aegilops_tauschii       AET0Gv20114300  15.3846 47.1698 0
...


will be splitted and
and dumped to separate files, each representing a rice -> otherSpecies projection

=cut

my $file = shift;

my $p = new Text::RecordParser::Tab;
$p->comment( qr/^#/ ); 

$p->filename($file);
#$p->bind_header;
$p->bind_fields( qw[ species stable_id  otherspecies  other_stable_id   perc_id other_perc_id is_high_confidence] );
  
my $records = $p->fetchall_arrayref( { Columns => {} } );
  
my @records_array = @{$records};

print "Total number of records is ", scalar @records_array, "\n";

my %species_pair_hash;

shift @records_array;
foreach ( @records_array){

#print "spliced ", scalar @pair, "records\n";

    my $otherspecies = $_->{otherspecies};
    my $species = $_->{species};
    die ("ERROR: bad otherspecies $otherspecies, rice spcies is $species") if $otherspecies eq 'oryza_sativa';
    
    my $filename_base = camelize($otherspecies).'_osj';

    push @{$species_pair_hash{ $filename_base }}, join ("\t", ($_->{stable_id}, $_->{other_stable_id}, $_->{perc_id}, $_->{other_perc_id}, $_->{is_high_confidence}));
    
}

	
for my $filename_base( keys %species_pair_hash){

	open my $fh, '>', "${filename_base}.rtm" or die "cannot open file $filename_base.rmt to write";	 
	map{ print $fh "$_\n" } @{$species_pair_hash{$filename_base}};
	
	close $fh;
}

__END__

=head1 AUTHOR

   Sharon Wei (weix@cshl.edu)
   Gramene Project (www.gramene.org)
   Ware Lab
   Cold Spring Harbor Laboratory

   This is free software and may be distributed on the same terms as Perl itself.

=cut

