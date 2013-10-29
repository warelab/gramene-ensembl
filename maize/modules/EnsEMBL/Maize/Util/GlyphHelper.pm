package EnsEMBL::Maize::Util::GlyphHelper;

use strict;

use Readonly;

sub new {
    my $class = shift;
    my $self = bless {}, $class;

    return $self;
}

=head2 add_to_legend

    Adds legend values

=cut
sub add_to_legend {
    my $self = shift;
    my ($glyphset, $key, $colour) = @_;
    
    my $track = $glyphset->check();
    
    $glyphset->{'config'}->{'legends'} ||= +{};
    my $legend = $glyphset->{'config'}->{'legends'}->{$track} ||= +{
        'features' => {
            'priority' => 1,
            'legend'   => [],
        }
    };
    push_unique($legend->{'features'}->{'legend'}, $key => $colour);
}

=head2 push_unique

    Pushes unique values onto an array

=cut
sub push_unique {
    my ($arrayref, $key, $value) = @_;
    for (my $i = 0; $i < scalar @$arrayref; $i += 2) {
        if ($arrayref->[$i] eq $key) {
            return;
        }
    }
    push @$arrayref, ($key => $value);
}


1;
