package Bio::EnsEMBL::GlyphSet::quality_score;

=head1 NAME

EnsEMBL::Web::GlyphSet::quality_score;

=head1 SYNOPSIS

Render wiggle plot of quality scores

=head1 LICENCE

This code is distributed under an Apache style licence:
Please see http://www.ensembl.org/code_licence.html for details

=head1 CONTACT

Shiran Pasternak - shiran@cshl.edu

=cut

use strict;
use vars qw(@ISA);
use Bio::EnsEMBL::GlyphSet;
use Sanger::Graphics::Glyph::Rect;
use Sanger::Graphics::Glyph::Composite;
use Sanger::Graphics::Bump;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::GlyphSet::wiggle_and_block;
use Time::HiRes qw(time);
use List::Util qw(max min);
use Memoize;

memoize('point_from_feature');

@ISA = qw(Bio::EnsEMBL::GlyphSet::wiggle_and_block);

sub init_label {
    my ($self) = @_;
    return if (defined $self->{'config'}->{'_no_label'});
    my $HELP_LINK = $self->check();
    $self->bumped(
        $self->{'config'}->get($HELP_LINK, 'compact') ? 'no' : 'yes')
        ;    # makes track expandable
    $self->init_label_text($self->my_config('label') || '---', $HELP_LINK);

    # $self->label->{zmenu} = { 'caption' => 'Options' };
    # if ($self->plottype eq 'area') {
    #     $self->label->{zmenu}->{'01:View as Line Plot'}
    #         = $self->url_with_param('opt_plottype', 'line');
    # } else {
    #     $self->label->{zmenu}->{'01:View as Area Plot'}
    #         = $self->url_with_param('opt_plottype', 'area');
    # }
    # if ($self->gridlines) {
    #     $self->label->{zmenu}->{'02:Hide Horizontal Gridlines'}
    #         = $self->url_with_param('opt_gridlines', 'off');
    # } else {
    #     $self->label->{zmenu}->{'02:Show Horizontal Gridlines'}
    #         = $self->url_with_param('opt_gridlines', 'on');
    # }
}

=head2 plottype

    Returns the plot type

=cut

sub plottype {
    my $self = shift;
    return $self->{config}->get('_settings', 'opt_plottype_line')
        ? 'line'
        : 'area';
}

=head2 gridlines

    Returns the gridlines setting

=cut

sub gridlines {
    my $self = shift;
    return $self->{config}->get('_settings', 'opt_gridlines');
}

=head2 url_with_param

    Returns a URL with the specified parameter value

=cut

sub url_with_param {
    my $self = shift;
    my ($param, $value) = @_;
    my $url = CGI::referer();
    $url =~ s/[;&]?$param=[^:&]*//;
    $url =~ s/[;&]$//;
    return "$url;$param=$value";
}

=head2 get_query_param

    Returns the parameter in the query

=cut

sub get_query_param {
    my $self    = shift;
    my ($param) = @_;
    my $query   = CGI->new(CGI::referer());
    return $query->param($param);
}

sub colour {
    return $_[0]->{'feature_colour'}, $_[0]->{'label_colour'},
        $_[0]->{'part_to_colour'};
}

sub draw_features {
    ### Called from {{ensembl-draw/modules/Bio/EnsEMBL/GlyphSet/wiggle_and_block.pm}}
    ### Arg 2 : draws wiggle plot if this is true
    ### Returns 0 if all goes well.
    ### Returns error message to print if there are features missing (string)

    my ($self, $db, $wiggle) = @_;
    my $type = $self->check();
    return unless defined $type;    ## No defined type arghhh!!

    my $strand      = $self->strand;
    my $Config      = $self->{'config'};
    my $strand_flag = $Config->get($type, 'str');
    return
        if ($strand_flag eq 'r' && $strand != -1
        || $strand_flag eq 'f' && $strand != 1);
    my $drawn_block = 0;
    my $container   = $self->{'container'};
    my $caption     = $Config->get($type, 'title')
        || $Config->get($type, 'label')
        || 'Quality scores';
    my %highlights;
    @highlights{ $self->highlights() } = ();
    my $length         = $container->length;
    my $pix_per_bp     = $Config->transform()->{'scalex'};
    my $DRAW_CIGAR     = $pix_per_bp > 0.2;
    my $feature_colour = $Config->get($type, 'col');
    my $h              = $Config->get('_settings', 'opt_halfheight') ? 4 : 8;
    my $chr            = $self->{'container'}->seq_region_name;

    my $X = -1e8;

    my @features = @{ $self->features($db) };

    $self->_offset($h);
    $self->render_track_name($caption, $feature_colour) if $drawn_block;

    my $drawn_wiggle = $wiggle ? $self->wiggle_plot($db) : 1;
    return 0 if $drawn_block && $drawn_wiggle;

    # Work out error message if some data is missing
    my $error;
    my $track = $self->{'config'}->get($type, 'label');

    if (!$drawn_block) {
        my $block_name = $self->my_config('block_name')
            || $self->my_config('label');
        $error .= $track;
    }

    if ($wiggle && !$drawn_wiggle) {
        $error .= " and " . $self->my_config('wiggle_name');
    }
    return $error;
}

sub features {
    ### Retrieves block features for constrained elements
    ### Returns arrayref of features

    my $self = shift;
    my ($db) = @_;

    my $slice                 = $self->{'container'};
    my $quality_score_adaptor = $db->get_adaptor('QualityScore');
    my $blocks = [];    #$quality_score_adaptor->fetch_all_by_Slice($slice);

    return $blocks;
}

sub wiggle_plot {

    ### Wiggle_plot
    ### Description: gets features for wiggle plot and passes to render_wiggle_plot
    ### Returns 1 if draws wiggles. Returns 0 if no wiggles drawn

    my ($self, $db) = @_;
    return 0 unless $db;

    $self->render_space_glyph();
    my $pix_per_base  = $self->{'config'}->transform()->{'scalex'};
    my $display_size  = $self->{'config'}->get('_settings', 'width') || 700;
    my $colour        = $self->my_config('col');
    my $display_label = 'Phred value';

    my $slice          = $self->{'container'};
    my $wiggle_adaptor = $db->get_QualityScoreAdaptor();
    if (!$wiggle_adaptor) {
        warn("Cannot get get adaptors: $wiggle_adaptor");
        return 0;
    }
    my $features
        = $wiggle_adaptor->fetch_all_by_Slice($slice, undef, $display_size)
        || [];
    return 0 unless scalar @$features;

    for my $feature (@$features) {
        $feature->start($feature->position - $slice->start + 1);
        $feature->end($feature->start + 1);
    }

    my $tick_marks = [ 0, 30, 60, 90 ];
    $self->render_wiggle_plot($features, $display_label, $tick_marks);

    return 1;
}

sub render_wiggle_plot {
    my ($self, $features, $display_label, $tick_marks) = @_;
    my @intervals       = sort { $a <=> $b } @$tick_marks;
    my $min_score       = $intervals[0];
    my $max_score       = $intervals[$#intervals];
    my $row_height      = 60;
    my $P_MAX           = max($max_score, 0);
    my $N_MIN           = min($min_score, 0);
    my $pix_per_score   = $row_height / ($P_MAX - $N_MIN);
    my $y_scale_factor  = $P_MAX / 100.0;
    my $red_line_offset = $P_MAX * $pix_per_score;
    my $pix_per_bp      = $self->{config}->transform->{'scalex'};
    my $tick_mark_width = 4;
    my ($fontname, $fontsize) = $self->get_font_details('innertext');

    my $SETTINGS = {
        'max_score'        => $P_MAX,
        'min_score'        => $N_MIN,
        'pixels_per_score' => $pix_per_score,
        'pixels_per_base'  => $pix_per_bp,
        'colour'           => $self->my_config('col'),
        'axis_colour'      => $self->my_config('axis_colour'),
        'grid_colour'      => $self->my_config('grid_colour'),
        'tick_mark_width'  => 4,
        'row_height'       => $row_height,
        'intervals'        => \@intervals,
        'y_factor'         => $y_scale_factor,
        'slice_length'     => $self->{container}->length(),
        'x_axis_offset'    => $red_line_offset,
        'inner_fontname'   => $fontname,
        'inner_fontsize'   => $fontsize,
        'axis_title'       => $display_label,
    };

    $self->render_tick_marks($SETTINGS);
    $self->render_quality_features($features, $SETTINGS);
    $self->render_axes($SETTINGS);

    return 1;
}

=head2 render_quality_features

    Add quality features to render

=cut

sub render_quality_features {
    my $self = shift;
    my ($features, $SETTINGS) = @_;

    my $get_height = undef;
    my $post_feature_render;

    if ($self->plottype eq 'line') {
        $get_height = sub { return 1; };
        $post_feature_render = sub {
            my ($point, $features, $i) = @_;
            return if ($i >= scalar @$features - 1);
            my $next_point
                = $self->point_from_feature($features->[ $i + 1 ], $SETTINGS);
            if ($point->{score} > 0 && $next_point->{score} > 0) {
                $self->push(
                    new Sanger::Graphics::Glyph::Line(
                        {   'pixelx'     => $point->{x} + $point->{width},
                            'pixelwidth' => $next_point->{x} 
                                - $point->{x}
                                - $point->{width},
                            'y'           => $point->{y},
                            'pixelheight' => $next_point->{y} - $point->{y},
                            'colour'      => $SETTINGS->{colour},
                        }
                    )
                );
            }
        };
    } else {
        $get_height          = sub { return $_[0]->{height} };
        $post_feature_render = sub { };
    }

    for (my $i = 0; $i < scalar @$features; $i++) {
        my $feature = $features->[$i];
        my $point = $self->point_from_feature($feature, $SETTINGS);
        $self->push(
            new Sanger::Graphics::Glyph::Rect(
                {   'y'         => $point->{y},
                    'height'    => $get_height->($point),
                    'x'         => $point->{x},
                    'width'     => $point->{width},
                    'absolutey' => 1,
                    'title'     => $point->{score},
                    'colour'    => $SETTINGS->{colour},
                    'zmenu'     => { 'caption' => "Phred $point->{score}" },
                }
            )
        );
        $post_feature_render->($point, $features, $i);
    }
}

=head2 point_from_feature

    Creates a drawable point from a score feature

=cut

sub point_from_feature {
    my $self = shift;
    my ($feature, $SETTINGS) = @_;
    my $point = +{};
    $point->{score} = $feature->score || 0;
    $point->{x} = max($feature->start, 0);
    $point->{width}
        = min($SETTINGS->{slice_length}, $feature->end) - $point->{x};
    $point->{height} = abs($point->{score} * $SETTINGS->{pixels_per_score});
    $point->{y}      = (
        $point->{score} < 0
        ? 0
        : -$point->{score} * $SETTINGS->{pixels_per_score}
    ) + $self->_offset + $SETTINGS->{x_axis_offset};
    return $point;
}

=head2 render_axes

    Render the axes

=cut

sub render_axes {
    my $self = shift;
    my ($SETTINGS) = @_;

    # X AXIS
    $self->push(
        new Sanger::Graphics::Glyph::Line(
            {   'x'         => 0,
                'y'         => $self->_offset + $SETTINGS->{x_axis_offset},
                'width'     => $SETTINGS->{slice_length},
                'height'    => 0,
                'absolutey' => 1,
                'colour'    => $SETTINGS->{axis_colour},
                'dotted'    => 1,
            }
        )
    );

    # Y AXIS
    $self->push(
        new Sanger::Graphics::Glyph::Line(
            {   'x'         => 0,
                'y'         => $self->_offset,
                'width'     => 0,
                'height'    => $SETTINGS->{row_height},
                'absolutey' => 1,
                'absolutex' => 1,
                'colour'    => $SETTINGS->{axis_colour},
                'dotted'    => 1,
            }
        )
    );
    $self->_offset($SETTINGS->{row_height});
    $self->render_x_axis_title($SETTINGS);
}

=head2 render_tick_marks

    Renders tick marks

=cut

sub render_tick_marks {
    my $self = shift;
    my ($SETTINGS) = @_;

    # Draw score tick marks
    for my $interval (@{ $SETTINGS->{intervals} }) {
        my $y
            = $self->_offset
            + ($SETTINGS->{max_score} - $interval)
            * $SETTINGS->{pixels_per_score};
        my @res_i = $self->get_text_width(
            0, $interval, '',
            'font'   => $SETTINGS->{inner_fontname},
            'ptsize' => $SETTINGS->{inner_fontsize}
        );
        my $textheight_i = $res_i[3];
        my $center_line  = $textheight_i / 2;
        $self->push(
            new Sanger::Graphics::Glyph::Text(
                {   'text'      => $interval,
                    'width'     => $res_i[2],
                    'textwidth' => $res_i[2],
                    'font'      => $SETTINGS->{inner_fontname},
                    'ptsize'    => $SETTINGS->{inner_fontsize},
                    'halign'    => 'right',
                    'valign'    => 'center',
                    'colour'    => $SETTINGS->{axis_colour},
                    'height'    => $textheight_i,
                    'y'         => $y - $center_line,
                    'x'         => -$SETTINGS->{tick_mark_width} - $res_i[2],
                    'absolutey' => 1,
                    'absolutex' => 1,
                    'absolutewidth' => 1,
                }
            )
        );

        # Draw the tick
        $self->push(
            new Sanger::Graphics::Glyph::Line(
                {   'x'             => -$SETTINGS->{tick_mark_width},
                    'y'             => $y,
                    'width'         => $SETTINGS->{tick_mark_width},
                    'height'        => 0,
                    'absolutey'     => 1,
                    'absolutex'     => 1,
                    'absolutewidth' => 1,
                    'colour'        => $SETTINGS->{axis_colour},
                }
            )
        );

        # Avoid rendering gridline at X axis
        if (   $self->gridlines
            && $y != $self->_offset + $SETTINGS->{x_axis_offset})
        {
            $self->push(
                new Sanger::Graphics::Glyph::Line(
                    {   'x'         => 0,
                        'y'         => $y,
                        'width'     => $SETTINGS->{slice_length},
                        'height'    => 0,
                        'absolutey' => 1,
                        'absolutex' => 1,
                        'colour'    => $SETTINGS->{grid_colour},
                    }
                )
            );
        }
    }
}

=head2 render_x_axis_title

    Adds a title to the X axis

=cut

sub render_x_axis_title {
    my $self         = shift;
    my ($SETTINGS)   = @_;
    my @res_analysis = $self->get_text_width(
        0, $SETTINGS->{axis_title},
        '',
        'font'   => $SETTINGS->{inner_fontname},
        'ptsize' => $SETTINGS->{inner_fontsize}
    );

    $self->push(
        new Sanger::Graphics::Glyph::Text(
            {   'text'      => $SETTINGS->{axis_title},
                'width'     => $res_analysis[2],
                'font'      => $SETTINGS->{inner_fontname},
                'ptsize'    => $SETTINGS->{inner_fontsize},
                'halign'    => 'left',
                'valign'    => 'bottom',
                'colour'    => $SETTINGS->{colour},
                'y'         => $self->_offset,
                'height'    => $res_analysis[3],
                'x'         => 1,
                'absolutey' => 1,
                'absolutex' => 1,
            }
        )
    );
    $self->_offset($res_analysis[3]);
}

1;
