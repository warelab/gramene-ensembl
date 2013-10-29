package EnsEMBL::Maize::Configuration::Chromosome;
use strict;
use EnsEMBL::Web::Configuration;
use EnsEMBL::Maize::Configuration;
our @ISA = qw( EnsEMBL::Web::Configuration );

## Function to configure map view
## MapView uses a single panel to display a chromosome image plus
## table rows containing basic information and 'navigation' forms

sub mapview {
    my $self = shift;

    my $document = $self->{page};
    my $content  = $document->content;
    my $panel    = $content->panel('info');
    
    # make cacheable --> didn't work
#    $panel->{'cacheable'} = 'yes';
    
    $panel->remove_component('change_chr');    # Move to sidebar

    # remove contig view to make way for fpc ctg and accessioned bac pulldowns
    $panel->remove_component('jump_to_contig');

    # changed method from add_component_first
    $panel->add_component(
        'disclaimer' => 'EnsEMBL::Maize::Component::Chromosome::disclaimer');

    $panel->replace_component(
        'image' => 'EnsEMBL::Maize::Component::Chromosome::chr_map');

    $panel->replace_component(
        'stats' => 'EnsEMBL::Maize::Component::Chromosome::stats');
    
    $panel->add_component('jump_to_contig' =>
            'EnsEMBL::Maize::Component::GenomeEntryPoints::fpc_form');

    $panel->add_component(
        'jump_to_bac' => 'EnsEMBL::Maize::Component::GenomeEntryPoints::bac_form'
    );

    my $chromosome = $self->{object};
    my $species    = $chromosome->species;
    my $chr_name   = $chromosome->chr_name;
    my @chr_options;

    # In maize context, the switch chrom context menu option is redundant 
    # because of the navigation panel. Kept for reference.

#    foreach my $chr (@{ $chromosome->species_defs->ENSEMBL_CHROMOSOMES || [] })
#    {
#        push(
#            @chr_options,
#            {   text => "Chromosome $chr",
#                href => "/$species/mapview?chr=$chr",
#                raw  => 1
#            }
#        );
#    }

#    my $menu = $document->menu;
#    $menu->add_entry_after(
#        'chromosome', 'syntview',
#        code    => 'switch',
#        text    => 'Switch Chromosome',
#        href    => "/$species/mapview?chr=$chr_name",
#        options => [@chr_options]
#    );
}

# maize specific mapview context menu 

sub context_menu {
    my $self = shift;
    my $obj      = $self->{object};
    my $species  = $obj->species;
    my $chr_name = $obj->chr_name;
    
    my $flag     = "chromosome";
    
    # remove menu redundancy with really_delete_menu_block
    my $menu = $self->{page}->menu;
    EnsEMBL::Maize::Configuration::really_delete_menu_block($menu, $flag);

    
    $self->{page}->menu->add_block( $flag, 'bulleted', "Chromosome $chr_name" );
  # create synteny form if relevant                                                                                                                                       
    my %hash  = $obj->species_defs->multi('SYNTENY');
    my @SPECIES = grep { @{ $obj->species_defs->other_species( $_, 'ENSEMBL_CHROMOSOMES' )||[]} } keys( %hash );
    
    # redirecting to cytoview 
    
    if( $chr_name ) {
	$self->{page}->menu->add_entry( $flag, code => 'name', 'text' => "@{[$obj->seq_region_type_and_name]} in cytoview",
					'href' => "/$species/cytoview?chr=$chr_name" );
    }
    if( @SPECIES ){
	my $array_ref;
	foreach my $next (@SPECIES) {
	    my $bio_name = $obj->species_defs->other_species($next, "SPECIES_BIO_NAME");
	    my $common_name = $obj->species_defs->other_species($next, "SPECIES_COMMON_NAME");
	    my $hash_ref = {'text'=>"vs $common_name (<i>$bio_name</i>)", 'href'=>"/$species/syntenyview?otherspecies=$next;chr=$chr_name", 'raw'=>1} ;
	    push (@$array_ref, $hash_ref);
	}

	$self->{page}->menu->add_entry($flag, code => 'syntview', 'href' => $array_ref->[0]->{'href'}, 'text'=>"View Chr $chr_name Synteny",
    'options' => $array_ref,
				       );

    }
    $self->{page}->menu->add_entry( $flag, code => 'karview', 'text' => "Map your data onto this chromosome",
				    'href' => "/$species/karyoview?chr=$chr_name" );
}


#----------------------------------------------------------------------
# Maize-specific syntenyview config. No need for homology yet.
# Also - move chr selection to sidebar
sub syntenyview {
    my $self = shift;

    my $document = $self->{page};
    my $content  = $document->content;
    my $panel    = $content->panel('info');

    # Remove unwanted components
    $panel->remove_component('syn_matches');     # No homology yet
    $panel->remove_component('nav_homology');    # No homology yet
    $panel->remove_component('change_chr');      # Move to sidebar

    # Chromosome-picker as part of sidebar
    my $chromosome   = $self->{object};
    my $otherspecies = $chromosome->param("otherspecies");
    my $species      = $chromosome->species;
    my $chr_name     = $chromosome->chr_name;
    my @chr_options;

    foreach my $chr (@{ $chromosome->species_defs->ENSEMBL_CHROMOSOMES || [] })
    {
        push(
            @chr_options,
            {   text => "Chromosome $chr",
                href => (
                          "/$species/syntenyview"
                        . "?chr=$chr&otherspecies=$otherspecies"
                ),
                raw => 1
            }
        );
    }

    my $menu  = $document->menu;
    my $block = $menu->block('chromosome');
    $menu->add_entry_after(
        'chromosome',
        'syntview',
        code => 'switch',
        text => 'Switch Chromosome',
        href => (
                  "/$species/syntenyview"
                . "?chr=$chr_name&otherspecies=$otherspecies"
        ),
        options => [@chr_options]
    );

}

1;
