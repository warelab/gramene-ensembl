package EnsEMBL::Maize::Component::GenomeEntryPoints;

use strict;
use warnings;

use EnsEMBL::Web::Form;
use EnsEMBL::Maize::Util::FPC;
use Data::Dumper qw(Dumper);    # For debug

our @ISA = qw( EnsEMBL::Web::Component );

sub bac_form {
    my $panel = shift;
    my ($object) = shift || "null";

    my $fpcutil = EnsEMBL::Maize::Util::FPC->new;
    my $bacs    = $fpcutil->fetch_analyzed_bacs($object);
    my $label   = 'View BAC in ContigView';
    if (defined $object && $object->can('chr_name')) {
        $label
            = "View analyzed BAC in ContigView (@{[scalar @$bacs]} available)";
    }
    my $form = _get_drop_down_form(
        {   'formname'    => 'bac',
            'action'      => "/Zea_mays2/contigview",
            'attribute'   => 'clone',
            'label'       => $label,
            'values'      => $bacs,
            'nameisvalue' => 0,
        }
    );
    $panel->add_row($form->render);
}

sub fpc_form {
    my $panel = shift;
    my ($object) = @_;

    my $fpcutil = EnsEMBL::Maize::Util::FPC->new;
    my $contigs = $fpcutil->fetch_contigs($object);
    my $label   = 'View FPC contig in CytoView';
    my $form    = _get_drop_down_form(
        {   'formname'    => 'fpc',
            'action'      => '/Zea_mays/cytoview',
            'attribute'   => 'contig',
            'label'       => $label,
            'values'      => $contigs,
            'nameisvalue' => 1,
        }
    );
    $panel->add_row($form->render);
}

sub virtualcorebin_form {
    my $panel = shift;
    my ($object) = @_;

    my $label = 'View virtual bin in CytoView';

    my $fpcutil     = EnsEMBL::Maize::Util::FPC->new;
    my $virtualbins = $fpcutil->fetch_corebins($object);

    my @virtualbins = sort { $a <=> $b } keys %$virtualbins;
    my $form = _get_drop_down_form(
        {   'formname'    => 'virtualcorebin',
            'action'      => '/Zea_mays/cytoview',
            'attribute'   => 'mapfrag',
            'label'       => $label,
            'values'      => \@virtualbins,
            'nameisvalue' => 1,
        }
    );
    $panel->add_row($form->render);
}

sub corebinmarker_form {
    my $panel = shift;
    my ($object) = @_;

    my $fpcutil        = EnsEMBL::Maize::Util::FPC->new;
    my $corebinmarkers = $fpcutil->fetch_corebinmarkers($object);

    my $label = 'View core marker in CytoView';

    my $form = _get_drop_down_form(
        {   'formname'    => 'corebinmarker',
            'action'      => '/Zea_mays/cytoview',
            'attribute'   => 'marker',
            'label'       => $label,
            'values'      => $corebinmarkers,
            'nameisvalue' => 1,
        }
    );
    $panel->add_row($form->render);
}

sub mapview_form {
    my $panel = shift;
    my ($object) = @_;

    my $text = <<HTML;
    <div class="forminline">
<form>
<h6>MapView</h6>
<img src=
"/images/zmays_karyo.png" usemap="#zmays_karyo" style=
"border: none" alt="Maize Karyotypes" /> <map name="zmays_karyo"
id="zmays_karyo">
    <area shape="rect" coords="7,4,19,155" alt="Chromosome 1"
    title="Chromosome 1" href="/Zea_mays/mapview?chr=1" />
    <area shape="rect" coords="24,23,36,155" alt="Chromosome 2"
    title="Chromosome 2" href="/Zea_mays/mapview?chr=2" />
    <area shape="rect" coords="41,25,53,155" alt="Chromosome 3"
    title="Chromosome 3" href="/Zea_mays/mapview?chr=3" />
    <area shape="rect" coords="58,14,70,155" alt="Chromosome 4"
    title="Chromosome 4" href="/Zea_mays/mapview?chr=4" />
    <area shape="rect" coords="75,34,87,155" alt="Chromosome 5"
    title="Chromosome 5" href="/Zea_mays/mapview?chr=5" />
    <area shape="rect" coords="92,73,104,155" alt="Chromosome 6"
    title="Chromosome 6" href="/Zea_mays/mapview?chr=6" />
    <area shape="rect" coords="109,69,121,155" alt="Chromosome 7"
    title="Chromosome 7" href="/Zea_mays/mapview?chr=7" />
    <area shape="rect" coords="126,65,138,155" alt="Chromosome 8"
    title="Chromosome 8" href="/Zea_mays/mapview?chr=8" />
    <area shape="rect" coords="143,84,155,155" alt="Chromosome 9"
    title="Chromosome 9" href="/Zea_mays/mapview?chr=9" />
    <area shape="rect" coords="160,74,172,155" alt="Chromosome 10"
    title="Chromosome 10" href="/Zea_mays/mapview?chr=10" />
</map>
</form>
</div>
HTML
    $panel->add_row($text);
}

sub synteny_form {
    my $panel = shift;
    my ($object) = @_;
    my @chromosomes
        = map { +{ 'name' => "Chromosome $_", 'value' => $_ } } (1 .. 10);
    my $form = _get_drop_down_form(
        {   'formname'    => 'synteny',
            'action'      => '/Zea_mays/syntenyview',
            'attribute'   => 'chr',
            'label'       => 'View synteny with rice',
            'values'      => \@chromosomes,
            'nameisvalue' => 0,
            'hidden'      => [ [ 'otherspecies', 'Oryza_sativa' ], ],
        }
    );
    $panel->add_row($form->render);
}

sub _get_drop_down_form {
    my ($params) = @_;

    my @values =
        $params->{'nameisvalue'}
        ? map { +{ 'name' => $_, 'value' => $_ } } @{ $params->{'values'} }
        : @{ $params->{'values'} };

    my $form
        = EnsEMBL::Web::Form->new($params->{'formname'}, $params->{'action'},
        'get');
    $form->add_element(
        'select'   => 'select',
        'type'     => 'DropDown',
        'name'     => $params->{'attribute'},
        'label'    => $params->{'label'},
        'values'   => \@values,
        'spanning' => 'inline',
    );
    for my $hidden (@{ $params->{'hidden'} }) {
        $form->add_element(
            'type'  => 'Hidden',
            'name'  => $hidden->[0],
            'value' => $hidden->[1],
        );
    }

    $form->add_element(
        'type'     => 'Submit',
        'value'    => '>',
        'spanning' => 'inline',
    );
    return $form;
}

1;
