package EnsEMBL::Maize::Component::Location;

use strict;
use warnings;

use EnsEMBL::Web::Component;
our @ISA = qw( EnsEMBL::Web::Component);

sub contigviewbottom_menu {
    my ($panel, $object) = @_;
    my $mc = $object->new_menu_container(
        'configname' => 'contigviewbottom',
        'panel'      => 'bottom',
        'leftmenus'  => [
            (   'Features',
                'ESTFeatures',
                'GSSFeatures',
                'FSTFeatures',
                'ArrayFeatures',
                'Repeats',
                'MDRFeatures',
                'Markers',
                'Compara',
                'DAS',
                'Options',
                'ImageSize',
            )
        ],
        'rightmenus' => [qw(Help)],
    );
    $panel->print($mc->render_html);
    $panel->print($mc->render_js);
    return 0;
}

1;
