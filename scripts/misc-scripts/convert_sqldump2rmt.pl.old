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

mysql-cabot -q ensembl_compara_plants_41_94 -e "SELECT genome_db.name as species, stable_id, perc_id FROM method_link_species_set INNER JOIN species_set USING(species_set_id) INNER JOIN homology USING(method_link_species_set_id) INNER JOIN homology_member USING(homology_id) INNER JOIN gene_member USING(gene_member_id) INNER JOIN genome_db ON gene_member.genome_db_id = genome_db.genome_db_id WHERE method_link_id = 201 AND species_set.genome_db_id = 1985" >Osj2allOrthologs.txt 

and looks like

species stable_id       perc_id
arabidopsis_thaliana    AT4G21870       20.8955
oryza_sativa    Os10g0437700    17.1779
arabidopsis_thaliana    AT4G21870       41.791
oryza_sativa    Os07g0517100    32.3699
...


will be converted to 

Os10g0437700	AT4G21870 17.1779 20.8955 1
Os07g0517100	AT4G21870	32.3699	41.791	1
...

merge every two lines into one line
and dumped to separate files, each representing a rice -> otherSpecies projection

=cut

my $file = shift;

my $p = new Text::RecordParser::Tab;
$p->comment( qr/^#/ ); 

$p->filename($file);
$p->bind_fields( qw[ species stable_id       perc_id ] );
  
my $records = $p->fetchall_arrayref( { Columns => {} } );
  
my @records_array = @{$records};

print "Total number of records is ", scalar @records_array, "\n";

my %species_pair_hash;
my @pair;

while ( @pair = splice (@records_array, 0, 2) ){

map{ warn ($_->{'species'}, ', ', $_->{'stable_id'}) } @pair;
    
#print "spliced ", scalar @pair, "records\n";
    warn "ERROR: less than two records in pair\n" if (scalar @pair < 2);

    if (scalar @pair < 2){
	warn "ERROR: less than two records in pair\n";
	map{ warn ($_->{'species'}, ', ', $_->{'stable_id'}) } @pair;
	next;
    }

    my ($rice_record, $other_record)  = $pair[1]->{'species'} eq 'oryza_sativa' ?
	($pair[1], $pair[0]) : ($pair[0], $pair[1]);

    my $filename_base = camelize($other_record->{'species'}).'_osj';

    push @{$species_pair_hash{ $filename_base }}, [$rice_record->{stable_id}, $other_record->{stable_id}, $rice_record->{perc_id}, $other_record->{perc_id}];
    
}

	
for my $filename_base( keys %species_pair_hash){

	open my $fh, '>', "${filename_base}.rmt" or die "cannot open file $filename_base.rmt to write";	 
	map{ print $fh join ("\t", (@{$_}, "1\n")) } @{$species_pair_hash{$filename_base}};
	
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

