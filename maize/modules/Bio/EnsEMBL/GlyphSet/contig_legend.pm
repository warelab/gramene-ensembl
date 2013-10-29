package Bio::EnsEMBL::GlyphSet::contig_legend;

use strict;
use Bio::EnsEMBL::GlyphSet;
use Sanger::Graphics::Glyph::Rect;
use Sanger::Graphics::Glyph::Text;
use Sanger::Graphics::Glyph::Poly;
use Bio::EnsEMBL::Utils::Eprof qw(eprof_start eprof_end);
our @ISA = qw(Bio::EnsEMBL::GlyphSet);

sub init_label {
    my ($self) = @_;
    return if (defined $self->{'config'}->{'_no_label'});
    $self->init_label_text('Contig legend');
}

sub _init {
    my ($self) = @_;

    return unless ($self->strand() == -1);

    my $BOX_WIDTH       = 20;
    my $NO_OF_COLUMNS   = 2;
    my $SEPARATOR_WIDTH = 3;
    my $KEY_WIDTH       = 10;

    my $vc       = $self->{'container'};
    my $Config   = $self->{'config'};
    my $im_width = $Config->image_width();
    my $type     = $Config->get('contig_legend', 'src');

    my @colours;
    return unless $Config->{'contig_legend_features'};
    my %features = %{ $Config->{'contig_legend_features'} };
    return unless %features;

    # Set up a separating line...
    my $rect = new Sanger::Graphics::Glyph::Rect(
        {   'x'             => 0,
            'y'             => 0,
            'width'         => $im_width,
            'height'        => 0,
            'colour'        => 'grey50',
            'absolutey'     => 1,
            'absolutex'     => 1,
            'absolutewidth' => 1,
        }
    );
    $self->push($rect);

    my ($x, $y) = (0, 0);
    my ($fontname, $fontsize) = $self->get_font_details('legend');
    my @res = $self->get_text_width(
        0, 'X', '',
        'font'   => $fontname,
        'ptsize' => $fontsize
    );
    my $th           = $res[3];
    my $pix_per_bp   = $self->{'config'}->transform()->{'scalex'};
    my $legend_width = $im_width / $NO_OF_COLUMNS;

    foreach (
        sort { $features{$b}->{'priority'} <=> $features{$a}->{'priority'} }
        keys %features
        )
    {
        @colours = @{ $features{$_}->{'legend'} };
        $y++ unless $x == 0;
        $x = 0;
        while (my ($legend_text, $legend) = splice @colours, 0, 2) {
            my $box_x_offset = 0;
            my $top          = $y * ($th + 3) + 2;
            my $bottom       = $top + $th - 2;
            if ($legend->{'type'} eq 'clone_end') {
                $self->push(
                    Bio::EnsEMBL::GlyphSet::contig::clone_end_glyph({
                        'clone_end' => 'LEFT',
                        'vector'    => $legend->{'value'},
                        'x_center' => $x * $legend_width + $box_x_offset + 5,
                        'y_start' => $top + 1,
                        'height'  => 10,
                    })
                );
                $box_x_offset += 15;
            } elsif ($legend->{'type'} eq 'colours') {
                for my $colour_pair (@{ $legend->{'value'} }) {
                    my $left   = $x * $legend_width + $box_x_offset;
                    my $right  = $left + $KEY_WIDTH;
                    my @points = (
                        $left, $top,    $right, $top,
                        $left, $bottom, $right, $bottom
                    );
                    for (my $i; $i < scalar @$colour_pair; $i++) {
                        my @triangle = @points[ ($i * 2) .. ($i * 2 + 5) ];
                        $self->push(
                            new Sanger::Graphics::Glyph::Poly(
                                {   'points'    => [@triangle],
                                    'colour'    => $colour_pair->[$i],
                                    'absolutex' => 1,
                                    'absolutey' => 1
                                }
                            )
                        );
                    }

                    $box_x_offset += $KEY_WIDTH + $SEPARATOR_WIDTH;
                }
            }
            $self->push(
                new Sanger::Graphics::Glyph::Text(
                    {   'x'             => $legend_width * $x + $box_x_offset,
                        'y'             => $y * ($th + 3) + 1,
                        'height'        => $th,
                        'valign'        => 'center',
                        'halign'        => 'left',
                        'ptsize'        => $fontsize,
                        'font'          => $fontname,
                        'colour'        => 'black',
                        'text'          => $legend_text,
                        'absolutey'     => 1,
                        'absolutex'     => 1,
                        'absolutewidth' => 1,
                    }
                )
            );
            $x++;
            if ($x == $NO_OF_COLUMNS) {
                $x = 0;
                $y++;
            }
        }
    }
}

1;

