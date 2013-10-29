
package EnsEMBL::Maize::Document::HTML::GenomeEntryPoints;

use strict;
use warnings;

use EnsEMBL::Web::Document::Panel;

our @ISA = qw( EnsEMBL::Web::Document::Panel );

sub new {
    my $class    = shift;
    my $self     = $class->SUPER::new();
    my ($params) = @_;

    $self->renderer = new EnsEMBL::Web::Document::Renderer::String();

    if ($params->{'all'} == 1) {
        $self->add_all_components();
    }

    return $self;
}

sub add_all_components {
    my $self = shift;
    for my $form_type (qw(fpc virtualcorebin corebinmarker mapview synteny))
    {
        $self->add_component($form_type);
    }
}

sub add_fpc_component {
    my $self = shift;
    my ($params) = @_;
    $self->add_component('fpc');
}

sub add_component {
    my $self   = shift;
    my ($type) = @_;
    my $form   = "${type}_form";
    $self->SUPER::add_component(
        $form => "EnsEMBL::Maize::Component::GenomeEntryPoints::$form");
}

sub add_row {
    my $self = shift;
    my ($content) = @_;
    $self->print(
        qq{
<tr>
    <td>$content</td>
</tr>
    }
    );
}

sub _start {
    $_[0]->print(
        q{
<div id="navigation_div">
<table id="genome_navigation" border="0" cellpadding="0" cellspacing="0">
    }
    );
}

sub _end {
    $_[0]->print(
        q{
</table>
</div>
    }
    );
}

1;
