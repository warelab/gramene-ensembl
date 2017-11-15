#! /usr/bin/env perl

#cat tree_id.* | ./script.taxon_pie.pl > taxon_pie
# Create annotation data for iTOL pie charts, see "http://itol.embl.de/help/dataset_piechart_template.txt"
# Output will be copied into laptop file:/Users/steinj/Documents/GRAMENE/Release_52/tree_data/iTOL_annotation/taxon_pies.txt

# Input taxa:
#Andropogoneae
#BEP_clade
#Brachypodium distachyon
#Mesangiospermae
#N
#Oryza sativa Japonica Group
#Panicoideae
#Poaceae
#Setaria italica
#Sorghum bicolor
#zea
#Zea mays

use strict;

my (%taxa, %tax_cnt);
my %colors = get_colors();
my %sp_convert = species_convert();
my %tax_convert = taxon_convert();
my $species;

while(<>){
    chomp;
    next if /ref_name/;
    if (/species:/){
	$species = $_; 
	$species =~ s/#species://;
	$species = $sp_convert{$species} || $species;
	next;
    }
    my @f = split /\t/;
    my $tax = $f[7];
#    $tax =~ s/ /_/g;
    $tax = $tax_convert{$tax} || $tax;
    if (!$tax){
	$tax = 'genome-specific'; #$sp_convert{$species} || $species;
    }

    $taxa{$tax} = 1;
    $tax_cnt{$species}->{$tax}++;
#    print join ("\t", $f[0], $species, $tax),"\n";
}

my @taxa = (qw/Mesangiospermae Poaceae Panicoideae Andropogoneae BEP_clade Zea Zea_mays genome-specific/); #sort { $colors{$a} cmp $colors{$b} } keys %taxa;

print join(" ", 'FIELD_LABELS', @taxa), "\n"; 
print join(" ", 'FIELD_COLORS', map { $colors{$_} } @taxa, ), "\n";
print 'DATA', "\n";

for my $sp (sort keys %tax_cnt){
    my @line = ($sp, );
    for my $tax (@taxa){
	push @line, $tax_cnt{$sp}->{$tax} || 0;
    }
    print join(" ", @line),"\n";
}

sub taxon_convert {
    my %tax_convert = (
#	'BEP_clade' => 'BEP_clade',  #this is how it is written accidentially in the iTOL tree "NAM_grant";
	'Brachypodium distachyon' => 'genome-specific',
	'Oryza sativa Japonica Group' => 'genome-specific',
	'Setaria italica' => 'genome-specific',
	'Sorghum bicolor' => 'genome-specific',
	'zea_parviglumis' => 'genome-specific',
	'N' => 'Zea_mays',
	'zea' => 'Zea',
	'Zea mays' => 'genome-specific'
	);
    return %tax_convert;
}

sub species_convert {
    my %sp_convert = (
	'zea_mays' => 'Zea_mays_B73',
	'zea_maysph207' => 'Zea_mays_PH207',
	'zea_maysw22' => 'Zea_mays_W22',
	'zea_parviglumis' => 'Zea_parviglumis',
	'brachypodium_distachyon' => 'Brachypodium_distachyon',
	'oryza_sativa' => 'Oryza_sativa',
	'setaria_italica' => 'Setaria_italica',
	'sorghum_bicolor' => 'Sorghum_bicolor',
	);
    return %sp_convert;
}

sub get_colors {
    my %colors = (
	'Mesangiospermae' => '#aa58c3',
	'Poaceae' => '#4ca3e8',
	'Panicoideae' => '#f0a30d',
	'BEP_clade' => '#f0a30d',
	'Andropogoneae' => '#40bd38',
	'Zea' => '#675ad9',
	'Zea_mays' => '#b7c732',
	'genome-specific' => '#db3f3f',
	#'Zea_mays_W22' => '#db3f3f',
	#'Zea_mays_PH207' => '#db3f3f',
	#'Zea_parviglumis' => '#db3f3f',
	#'Brachypodium_distachyon' => '#db3f3f',
	#'Oryza_sativa' => '#db3f3f',
	#'Setaria_italica' => '#db3f3f',
	#'Sorghum_bicolor' => '#db3f3f',
	);
    return %colors;
}
