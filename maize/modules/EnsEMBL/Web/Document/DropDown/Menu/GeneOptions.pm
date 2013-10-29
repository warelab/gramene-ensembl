package EnsEMBL::Web::Document::DropDown::Menu::GeneOptions;

use strict;
use EnsEMBL::Web::Document::DropDown::Menu;

@EnsEMBL::Web::Document::DropDown::Menu::GeneOptions::ISA
    = qw( EnsEMBL::Web::Document::DropDown::Menu );

sub new {
    my $class = shift;
    my $self  = $class->SUPER::new(
        @_,    ## This contains the menu containers as the first element
        'image_name'  => 'y-gene_options',
        'image_width' => 57,
        'alt'         => 'Gene options'
    );

    my @menu_entries
        = @{ $self->{'config'}->get('_settings', 'gene_options') || [] };

    return undef unless @menu_entries;
    for my $menu_item (@menu_entries) {
        my ($key, $label) = @$menu_item;
        if ($self->{'config'}->is_setting($key)) {
            $self->add_checkbox(@$menu_item);
        }
    }
    return $self;
}

1;
