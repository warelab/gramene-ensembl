package Bio::EnsEMBL::GlyphSet::cyto_dna_alignments;

use strict;
use vars qw(@ISA);
use Bio::EnsEMBL::GlyphSet_simple;
@ISA = qw(Bio::EnsEMBL::GlyphSet_simple);

use EnsEMBL::Maize::Util::FPC;
use EnsEMBL::Maize::Util::GlyphHelper;
use Sanger::Graphics::Glyph::Poly;

sub my_label {
    my $self = shift;
    return $self->my_config('label');
}

sub features {
    my $self = shift;

    my $fpc  = EnsEMBL::Maize::Util::FPC->new;
    my $bacs = $fpc->fetch_sequenced_bacs_for_slice($self->{'container'});
    return $bacs;
}

## Link back to this page centred on the map fragment

sub href {
    my ($self, $f) = @_;
    return
        "/@{[$self->{container}{_config_file_name_}]}/$ENV{'ENSEMBL_SCRIPT'}?mapfrag=@{[$f->get_scalar_attribute('name')]}";
}

sub zmenu {
    my ($self, $feature, $sets, $features) = @_;
    my $subtype_labeler = $self->my_config('subtype_labeler');
    my $caption         = $self->my_label();
    my $zmenu = +{ 'caption' => "$caption: " . $features->{'All features'} };
    my $zmenu_index = 1;
    for my $subtype (@$sets) {
        if ($features->{$subtype}) {
            $zmenu->{ "$zmenu_index:"
                    . $subtype_labeler->($subtype)
                    . ": $features->{$subtype}" } = '';
            $zmenu_index++;
        }
    }
    return $zmenu;
}

=head2 colour

    The feature color

=cut

sub colour {
    my $self = shift;
    return $self->my_config('col');
}

=head2 feature_method

    Returns the name of the method to fetch features from a slice

=cut

sub feature_method {
    my $self = shift;
    return $self->my_config('feature_method') || 'get_all_DnaAlignFeatures';
}

sub _init {
    my ($self) = @_;

    use constant VSPACING => 2;

    my $fpc = EnsEMBL::Maize::Util::FPC->new();

    my $pixels_per_base = $self->{'config'}->transform()->{'scalex'};

    my $threshold    = $self->my_config('threshold');
    my $region_width = $self->{'container'}->length();

    if ($threshold < $region_width) {
        my $printable_threshold
            = $self->thousandify(sprintf('%d', $threshold));
        $self->errorTrack(
            sprintf(
                '%s only displayed for less than %s bp',
                $self->my_label, $printable_threshold
            )
        );
        return;
    }
    my $w            = $self->my_config('width');
    my $h            = $self->my_config('height');
    my $rw           = $w / 2;
    my $rh           = $h / 2;
    my $window_width = $region_width * $pixels_per_base;

    my $sets           = $self->my_config('sets');
    my $bac_features   = $self->features();
    my $argument_index = $self->my_config('argument_index') || 0;

    for my $feature (@$bac_features) {
        my $y = 0;
        my $alignments
            = $fpc->bac_dna_features($feature, $self->feature_method(), $sets,
            $argument_index);
        next
            unless (defined $alignments and scalar keys %$alignments > 0);

        my $colour = $self->colour();

        my $start = $feature->start;
        $start = 1 if ($start < 1);
        my $end = $feature->end;
        $end = $region_width if ($end > $region_width);
        my $c = ($end + $start) / 2 * $pixels_per_base;

        # Rhombus:
        my @polygon = (
            $c, $y,    # top center
            $c + $rw, $y + $rh,    # middle right
            $c, $y + $h,           # bottom center
            $c - $rw, $y + $rh     # middle left
        );

        # Triangle:
        # my @polygon = (
        #     $c, $y,    # top center
        #     $c + $rw, $y + $h,    # bottom right
        #     $c - $rw, $y + $h     # bottom left
        # );

        # Heuristic: avoid clipping pipeline if center point not near
        # the edge
        if ($c < $rw || $c > $window_width - $rw) {
            @polygon = $self->clip(\@polygon, $window_width);
        }

        $self->push(
            new Sanger::Graphics::Glyph::Poly(
                {   'points'    => [@polygon],
                    'colour'    => $colour,
                    'absolutex' => 1,
                    'absolutey' => 1,
                    'zmenu'     => $self->zmenu($feature, $sets, $alignments),
                    'href'      => $self->href($feature),
                }
            )
        );
        $y += $h + VSPACING;
    }
}

=head2 inside_clip_plane

    Determines if vertex is inside clip plane

=cut

sub inside_clip_plane {
    my $self = shift;
    my ($vertex, $length) = @_;

    return $vertex->[0] > 0 && $vertex->[0] < $length;
}

=head2 intersection

    Finds the vertex that intersects the line between two points in a given
    window

=cut

sub intersection {
    my $self = shift;
    my ($p1, $p2, $length) = @_;

    for my $x (0, $length) {
        if (   ($p1->[0] > $x && $p2->[0] < $x)
            || ($p1->[0] < $x && $p2->[0] > $x))
        {
            my $m = ($p2->[1] - $p1->[1]) / ($p2->[0] - $p1->[0]);
            my $y = $p2->[1] + $m * ($x - $p2->[0]);
            return [ $x, $y ];
        }
    }
    return [];
}

=head2 clip

    Clips a polygon if it intersects the window borders. Uses the
    Sutherland-Hodgman plane-clipping algorithm.

=cut

sub clip {
    my $self = shift;
    my ($polygon, $window_length) = @_;

    my @new_polygon = ();
    my @S           = @$polygon[ $#$polygon - 1 .. $#$polygon ];
    for (my $i = 0; $i < scalar @$polygon; $i += 2) {
        my @P = @$polygon[ $i, $i + 1 ];
        if ($self->inside_clip_plane(\@P, $window_length)) {
            if ($self->inside_clip_plane(\@S, $window_length)) {
                push @new_polygon, @P;
            } else {
                push @new_polygon,
                    @{ $self->intersection(\@S, \@P, $window_length) };
                push @new_polygon, @P;
            }
        } elsif ($self->inside_clip_plane(\@S, $window_length)) {
            push @new_polygon,
                @{ $self->intersection(\@P, \@S, $window_length) };
        }
        @S = @P;
    }
    return @new_polygon;
}

1;
