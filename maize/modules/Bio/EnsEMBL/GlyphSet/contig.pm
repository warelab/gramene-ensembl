package Bio::EnsEMBL::GlyphSet::contig;
use strict;
use vars qw(@ISA);
use Bio::EnsEMBL::GlyphSet_simple;
@ISA = qw(Bio::EnsEMBL::GlyphSet_simple);
use Readonly;
use Sanger::Graphics::Glyph::Rect;
use Sanger::Graphics::Glyph::Poly;
use Sanger::Graphics::Glyph::Space;
use Sanger::Graphics::Glyph::Text;
use Sanger::Graphics::ColourMap;
use Data::Dumper;
$Data::Dumper::Indent = 2;
use constant MAX_VIEWABLE_ASSEMBLY_SIZE => 5e6;

Readonly my $ODD_SCAFFOLD   => 1;
Readonly my $EVEN_SCAFFOLD  => 2;
Readonly my $NO_SCAFFOLD    => 3;
Readonly my %CONTIG_COLOURS => +{
    $ODD_SCAFFOLD  => [ 'contigblue1', 'contigblue2' ],
    $EVEN_SCAFFOLD => [ 'darkorange2', 'darkorange4' ],
    $NO_SCAFFOLD   => [ 'gray60',      'gray40' ],
};

Readonly my %VECTOR_COLOURS => (
    'SP6' => 'gold1',
    'T7'  => 'firebrick2',
);

my %color_cache = ();

sub init_label {
    my ($self) = @_;
    return if (defined $self->{'config'}->{'_no_label'});
    my $cs_name = 'contigs';
    if (my $tcs = $self->{container}->coord_system) {
        if (my $scs = $tcs->adaptor->fetch_sequence_level) {
            $cs_name = $scs->name;
        }
    }
    my ($fontname, $fontsize) = $self->get_font_details('label');
    my $label = new Sanger::Graphics::Glyph::Text(
        {   'text'      => "DNA:$cs_name",
            'font'      => $fontname,
            'ptsize'    => $fontsize,
            'absolutey' => 1,
            'href' =>
                "javascript:X=hw('@{[$self->{container}{_config_file_name_}]}','$ENV{ENSEMBL_SCRIPT}','contig')",
            'zmenu' => {
                'caption' => 'HELP',
                '01:Track information...' =>
                    "javascript:X=hw('@{[$self->{container}{_config_file_name_}]}','$ENV{ENSEMBL_SCRIPT}','contig')"
            }
        }
    );
    $self->label($label);
}

sub _init {
    my ($self) = @_;

    # only draw contigs once - on one strand
    return unless ($self->strand() == 1);

    my $Container = $self->{'container'};
    $self->{'vc'} = $Container;
    my $length = $Container->length();

    my $ystart = 3;

    my $gline = new Sanger::Graphics::Glyph::Rect(
        {   'x'         => 0,
            'y'         => $ystart + 7,
            'width'     => $length,
            'height'    => 0,
            'colour'    => 'grey50',
            'absolutey' => 1,
        }
    );

    $self->push($gline);

    my @features = ();
    foreach my $segment (@{ $Container->project('seqlevel') || [] }) {
        my $start     = $segment->from_start;
        my $end       = $segment->from_end;
        my $ctg_slice = $segment->to_Slice;
        my $ORI       = $ctg_slice->strand;
        my $feature   = {
            'start' => $start,
            'end'   => $end,
            'name'  => $ctg_slice->seq_region_name
        };
        $feature->{'locations'}{ $ctg_slice->coord_system->name } = [
            $ctg_slice->seq_region_name, $ctg_slice->start,
            $ctg_slice->end,             $ctg_slice->strand
        ];
        foreach (
            @{  $Container->adaptor->db->get_CoordSystemAdaptor->fetch_all()
                    || []
            }
            )
        {
            my $path;
            eval { $path = $ctg_slice->project($_->name); };
            next unless (@$path == 1);
            $path = $path->[0]->to_Slice;

            # get clone id out of seq_region_attrib for link to webFPC
            if ($_->{'name'} eq 'clone') {
                my ($clone_name)
                    = @{ $path->get_all_Attributes('fpc_clone_id') };
                $feature->{'internal_name'} = $clone_name->{'value'}
                    if $clone_name;
            }
            $feature->{'locations'}{ $_->name } = [
                $path->seq_region_name, $path->start,
                $path->end,             $path->strand
            ];
        }
        $feature->{'ori'} = $ORI;
        push @features, $feature;
    }

    if (@features) {
        $self->_init_non_assembled_contig($ystart, \@features);
    } else {
        $self->errorTrack("Golden path gap - no contigs to display!");
    }
}

sub _init_non_assembled_contig {
    my ($self, $ystart, $contig_tiling_path) = @_;

    my $Container = $self->{'vc'};

    #    warn Dumper ($Container);
    #    warn Dumper ($contig_tiling_path);
    my $length = $Container->length();
    my $ch     = $Container->seq_region_name;

    my $Config = $self->{'config'};

    my $module = ref($self);
    $module = $1 if $module =~ /::([^:]+)$/;
    my $threshold_navigation
        = ($Config->get($module, 'threshold_navigation') || 2e6) * 1001;
    my $navigation = $Config->get($module, 'navigation') || 'on';
    my $show_navigation
        = ($length < $threshold_navigation) && ($navigation eq 'on');

    my $cmap = $Config->colourmap();

    ########
    # Vars used only for scale drawing
    #
    my $black      = 'black';
    my $red        = 'red';
    my $highlights = join('|', $self->highlights());
    $highlights = $highlights ? "&highlight=$highlights" : '';
    if ($self->{'config'}->{'compara'})
    {    ## this is where we have to add in the other species....
        my $C = 0;
        foreach (@{ $self->{'config'}{'other_slices'} }) {
            if ($C != $self->{'config'}->{'slice_number'}) {
                if ($C) {
                    if ($_->{'location'}) {
                        $highlights .= sprintf(
                            "&s$C=%s&c$C=%s:%s:%s&w$C=%s",
                            $_->{'location'}->species,
                            $_->{'location'}->seq_region_name,
                            $_->{'location'}->centrepoint,
                            $_->{'ori'},
                            $_->{'location'}->length
                        );
                    } else {
                        $highlights .= sprintf("&s$C=%s", $_->{'species'});
                    }
                } else {
                    $highlights .= sprintf("&c=%s:%s:1&w=%s",
                        $_->{'location'}->seq_region_name,
                        $_->{'location'}->centrepoint,
                        $_->{'location'}->length);
                }
            }
            $C++;
        }
    }    ##

    my $contig_strand = $Container->can('strand') ? $Container->strand : 1;
    my $clone_based = $Config->get('_settings', 'clone_based') eq 'yes';
    my $global_start =
          $clone_based
        ? $Config->get('_settings', 'clone_start')
        : $Container->start();
    my $global_end = $global_start + $length - 1;
    my $im_width   = $Config->image_width();

    my $pixels_per_bp = $im_width / $length;

    #######
    # Draw the Contig Tiling Path
    #
    my ($w, $h) = $Config->texthelper()->real_px2bp('Tiny');
    my $i = 1;

    my $adaptor = $self->{'container'}->adaptor()->db()->get_SliceAdaptor();

    # Used for alternating colours
    my $colour_position = 0;
    foreach my $tile (sort { $a->{'start'} <=> $b->{'start'} }
        @{$contig_tiling_path})
    {
        my $rend   = $tile->{'end'};
        my $rstart = $tile->{'start'};
        my $rid    = $tile->{'name'};
        my $strand = $tile->{'ori'};
        $rstart = 1       if $rstart < 1;
        $rend   = $length if $rend > $length;

        # fetch the order and orientation attribs for this contig
        my $value;
        my $slice = $adaptor->fetch_by_region('contig', $rid);
        my $contig_color_index = $NO_SCAFFOLD;
        my $scaffold_name = unique_slice_attribute($slice, 'contig-scaffold')
            || 'NONE';
        my ($scaffold_number) = ($scaffold_name =~ m/Scaffold(\d+)/);
        if (defined $scaffold_number) {
            if ($scaffold_number % 2 == 0) {
                $contig_color_index = $EVEN_SCAFFOLD;
            } else {
                $contig_color_index = $ODD_SCAFFOLD;
            }
        }
        my $clone_end    = unique_slice_attribute($slice, 'clone-end');
        my $clone_vector = unique_slice_attribute($slice, 'clone-vector');

        my %contig_feature = (
            'x'     => $rstart - 1,
            'y'     => $ystart + 2,
            'width' => $rend - $rstart + 1,
        );
        push @{ $self->{'rendered_contigs'} }, \%contig_feature;
        my $colour = $color_cache{ $slice->get_seq_region_id }
            ||= $CONTIG_COLOURS{$contig_color_index}[$colour_position];
        my $glyph = create_contig_glyph(\%contig_feature, $colour);

        $colour_position = ($colour_position + 1) % 2;

        my $script  = $ENV{'ENSEMBL_SCRIPT'};
        my $caption = 'Centre on';
        if ($script eq 'multicontigview') {
            $script  = 'contigview';
            $caption = 'Jump to contigview';
        }
        if ($navigation eq 'on') {
            foreach (qw(chunk supercontig clone scaffold contig)) {
                if (my $Q = $tile->{'locations'}->{$_}) {
                    $glyph->{'href'}
                        = qq(/@{[$self->{container}{_config_file_name_}]}/$script?ch=$ch&region=$Q->[0]);
                }
            }
        }
        my $label = '';

        #	warn "$label\n";
        if ($show_navigation) {
            $glyph->{'zmenu'} = { 'caption' => $rid, };
            my $POS = 10;
            for (qw( contig clone supercontig scaffold chunk)) {
                if (my $Q = $tile->{'locations'}->{$_}) {
                    my $name = $Q->[0];
                    $name =~ s/\.\d+$// if $_ eq 'clone';
                    $label ||= $tile->{'locations'}->{$_}->[0];

                    # (Shiran) Kludge to parse out BAC name and 'Contig'
                    # prefix contig label so as to display only contig name
                    $label =~ s/^[A-Z]+\d+\.\d+\-Contig(\d+)$/$1/;

                    #                    warn "$label\n";
                    (my $T = ucfirst($_)) =~ s/contig/Contig/g;
                    $glyph->{'zmenu'}{"$POS:$T $name"} = ''
                        unless $_ eq 'contig';
                    $POS++;

                    #add links to Ensembl and FPC (vega danio)
                    if (/clone/) {
                        if (my $fpc_url = $self->ID_URL('FPC_ENSEMBL', $name))
                        {

                            # For gramene maize FPC link
                            $glyph->{'zmenu'}{"$POS:Jump to FPC browser"}
                                = $fpc_url;
                            $POS++;
                        }
                        my $ens_URL = $self->ID_URL('EGB_ENSEMBL', $name);
                        $glyph->{'zmenu'}{"$POS:View in Ensembl"} = $ens_URL
                            if $ens_URL;
                        $POS++;
                        my $internal_clone_name = $tile->{'internal_name'};
                        my $fpc_URL
                            = $self->ID_URL('FPC', $internal_clone_name);
                        $glyph->{'zmenu'}{"$POS:View in WebFPC"} = $fpc_URL
                            if $fpc_URL && $internal_clone_name;
                        $POS++;
                    }
                    $glyph->{'zmenu'}{"$POS:EMBL source file"}
                        = $self->ID_URL('EMBL', $name)
                        if /clone/;
                    $POS++;
                    $glyph->{'zmenu'}{"$POS:$caption $T"}
                        = qq(/@{[$self->{container}{_config_file_name_}]}/$script?ch=$ch&region=$name);
                    $POS++;
                    $glyph->{'zmenu'}{"$POS:Export this $T"}
                        = qq(/@{[$self->{container}{_config_file_name_}]}/exportview?tab=fasta&type=feature&ftype=$_&id=$name);
                    $POS++;
                    if (/contig/) {
                        $glyph->{'zmenu'}{"$POS:Scaffold: $scaffold_name"}
                            = '';
                        $POS++;

                        if (defined $clone_vector) {
                            $glyph->{'zmenu'}
                                {"$POS:Vector Side: $clone_vector"} = '';
                            $POS++;
                        }
                        if (defined $clone_end) {
                            $glyph->{'zmenu'}{"$POS:Clone End: $clone_end"}
                                = '';
                            $POS++;
                        }
                    }
                }
            }
        }
        $self->push($glyph);

        my $y_offset = 11;
        my $x_center = int(($rend + $rstart) / 2 * $pixels_per_bp);
        if (defined $clone_end) {
            $self->push(
                clone_end_glyph(
                    {   'vector'      => $clone_vector,
                        'clone_end'   => $clone_end,
                        'x_center'    => $x_center,
                        'y_start'     => $ystart + 2,
                        'height'      => $y_offset - 1,
                        'pixel_width' => ($rend - $rstart + 2)
                            * $pixels_per_bp,
                    }
                )
            );
        }

        # $label = $strand > 0 ? "$label >" : "< $label";
        #	warn "$label\n";
        my $text_colour = 'white';
        my ($fontname, $fontsize) = $self->get_font_details('innertext');

        my @label_dimension = $self->get_text_width(
            0, $label, $label,
            'font'   => $fontname,
            'ptsize' => $fontsize
        );
        if ($label_dimension[0]) {
            my $tglyph = new Sanger::Graphics::Glyph::Text(
                {   'x'         => $x_center - int($label_dimension[2] / 2),
                    'y'         => $ystart + $y_offset,
                    'font'      => $fontname,
                    'ptsize'    => $fontsize,
                    'colour'    => $text_colour,
                    'text'      => $label,
                    'absolutey' => 1,
                    'absolutex' => 1,
                }
            );
            $self->push($tglyph);
        }
    }

    ######
    # Draw the scale, ticks, red box etc
    #
    my $gline = new Sanger::Graphics::Glyph::Rect(
        {   'x'             => 0,
            'y'             => $ystart,
            'width'         => $im_width,
            'height'        => 0,
            'colour'        => $black,
            'absolutey'     => 1,
            'absolutex'     => 1,
            'absolutewidth' => 1,
        }
    );
    $self->unshift($gline);
    $gline = new Sanger::Graphics::Glyph::Rect(
        {   'x'             => 0,
            'y'             => $ystart + 15,
            'width'         => $im_width,
            'height'        => 0,
            'colour'        => $black,
            'absolutey'     => 1,
            'absolutex'     => 1,
            'absolutewidth' => 1,
        }
    );
    $self->unshift($gline);

    ## pull in our subclassed methods if necessary
    if ($self->can('add_arrows')) {
        $self->add_arrows($im_width, $black, $ystart);
    }

    my $tick;
    my $interval = int($im_width / 10);

    # correct for $im_width's that are not multiples of 10
    my $corr = int(($im_width % 10) / 2);
    for (my $i = 1; $i <= 9; $i++) {
        my $pos = $i * $interval + $corr;
        $self->unshift(
            new Sanger::Graphics::Glyph::Rect(
                {    # the forward strand ticks
                    'x'             => 0 + $pos,
                    'y'             => $ystart - 4,
                    'width'         => 0,
                    'height'        => 3,
                    'colour'        => $black,
                    'absolutey'     => 1,
                    'absolutex'     => 1,
                    'absolutewidth' => 1,
                }
            )
        );
        $self->unshift(
            new Sanger::Graphics::Glyph::Rect(
                {    # the reverse strand ticks
                    'x'             => 0 + $pos,
                    'y'             => $ystart + 16,
                    'width'         => 0,
                    'height'        => 3,
                    'colour'        => $black,
                    'absolutey'     => 1,
                    'absolutex'     => 1,
                    'absolutewidth' => 1,
                }
            )
        );
    }

    # The end ticks
    $self->unshift(
        new Sanger::Graphics::Glyph::Rect(
            {   'x'             => 0 + $corr,
                'y'             => $ystart - 2,
                'width'         => 0,
                'height'        => 1,
                'colour'        => $black,
                'absolutey'     => 1,
                'absolutex'     => 1,
                'absolutewidth' => 1,
            }
        )
    );

    # the reverse strand ticks
    $self->unshift(
        new Sanger::Graphics::Glyph::Rect(
            {   'x'             => $im_width - 1 - $corr,
                'y'             => $ystart + 16,
                'width'         => 0,
                'height'        => 1,
                'colour'        => $black,
                'absolutey'     => 1,
                'absolutex'     => 1,
                'absolutewidth' => 1,
            }
        )
    );

    my $Container_size_limit = $Config->get('_settings', 'default_vc_size');

    # only draw a red box if we are in contigview top and there is a
    # detailed display
    my $rbs = $Config->get('_settings', 'red_box_start');
    my $rbe = $Config->get('_settings', 'red_box_end');
    if ($Config->get('_settings', 'draw_red_box') eq 'yes') {

        # only draw focus box on the correct display...
        $self->unshift(
            new Sanger::Graphics::Glyph::Rect(
                {   'x'            => $rbs - $global_start,
                    'y'            => $ystart - 4,
                    'width'        => $rbe - $rbs + 1,
                    'height'       => 23,
                    'bordercolour' => $red,
                    'absolutey'    => 1,
                }
            )
        );
        $self->unshift(
            new Sanger::Graphics::Glyph::Rect(
                {   'x'            => $rbs - $global_start,
                    'y'            => $ystart - 3,
                    'width'        => $rbe - $rbs + 1,
                    'height'       => 21,
                    'bordercolour' => $red,
                    'absolutey'    => 1,
                }
            )
        );
    }

    my $width = $interval * ($length / $im_width);
    my $interval_middle = $width / 2;

    if ($navigation eq 'on') {
        foreach my $i (0 .. 9) {
            my $pos = $i * $interval;

            # the forward strand ticks
            $self->unshift(
                new Sanger::Graphics::Glyph::Space(
                    {   'x'             => 0 + $pos,
                        'y'             => $ystart - 4,
                        'width'         => $interval,
                        'height'        => 3,
                        'absolutey'     => 1,
                        'absolutex'     => 1,
                        'absolutewidth' => 1,
                        'href'          => $self->zoom_URL(
                            $Container->seq_region_name,
                            $interval_middle + $global_start,
                            $length,
                            1,
                            $highlights,
                            $self->{'config'}->{'slice_number'},
                            $contig_strand
                        ),
                        'zmenu' => $self->zoom_zmenu(
                            $Container->seq_region_name,
                            $interval_middle + $global_start,
                            $length,
                            $highlights,
                            $self->{'config'}->{'slice_number'},
                            $contig_strand
                        ),
                    }
                )
            );

            # the reverse strand ticks
            $self->unshift(
                new Sanger::Graphics::Glyph::Space(
                    {   'x'             => $im_width - $pos - $interval,
                        'y'             => $ystart + 16,
                        'width'         => $interval,
                        'height'        => 3,
                        'absolutey'     => 1,
                        'absolutex'     => 1,
                        'absolutewidth' => 1,
                        'href'          => $self->zoom_URL(
                            $Container->seq_region_name,
                            $global_end + 1 - $interval_middle,
                            $length,
                            1,
                            $highlights,
                            $self->{'config'}->{'slice_number'},
                            $contig_strand
                        ),
                        'zmenu' => $self->zoom_zmenu(
                            $Container->seq_region_name,
                            $global_end + 1 - $interval_middle,
                            $length,
                            $highlights,
                            $self->{'config'}->{'slice_number'},
                            $contig_strand
                        ),
                    }
                )
            );
            $interval_middle += $width;
        }
    }
    $self->add_legend();
}

sub add_legend {
    my $self = shift;
    $self->{'config'}->{'contig_legend_features'}->{'contigs'} = +{
        'priority' => 1000,
        'legend'   => [
            'Scaffolds of ordered contigs' => +{
                'type'  => 'colours',
                'value' => [
                    $CONTIG_COLOURS{$ODD_SCAFFOLD},
                    $CONTIG_COLOURS{$EVEN_SCAFFOLD},
                ]
            },
            'SP6 Vector' => +{
                'type'  => 'clone_end',
                'value' => 'SP6',
            },
            'Unordered contigs' => +{
                'type' => 'colours',
                'value' =>
                    [ $CONTIG_COLOURS{$NO_SCAFFOLD}, [ 'white', 'white' ] ],
            },
            'T7 Vector' => +{
                'type'  => 'clone_end',
                'value' => 'T7',
            },
        ],
    };
}

sub create_contig_glyph {
    my ($feature, $colour) = @_;
    my $glyph = new Sanger::Graphics::Glyph::Rect(
        {   %$feature,
            'height'    => 11,
            'colour'    => $colour,
            'absolutey' => 1,
        }
    );
    return $glyph;
}

sub clone_end_glyph {
    my ($params)    = @_;
    my $vector      = $params->{'vector'};
    my $clone_end   = $params->{'clone_end'};
    my $x_center    = $params->{'x_center'};
    my $y_start    = $params->{'y_start'};
    my $pixel_width = $params->{'pixel_width'};
    my $height      = $params->{'height'};

    my $half_width  = int(($pixel_width + 1) / 2);
    my $side_length = int(($height + 1) / 2);
    my $absolutex   = 1;
    my $absolutey   = 1;
    my $colour = $VECTOR_COLOURS{$vector} || undef;

    my $half_contig = 0;
    my $arrow_width = $side_length;
    if (defined $clone_end) {
        if ($clone_end eq 'LEFT') {
            $half_contig = -$half_width;
        } elsif ($clone_end eq 'RIGHT') {
            $half_contig = $half_width;
            $arrow_width = -$arrow_width;
        }
    }
    my $x_start      = $x_center + $half_contig;
    my $x_arrow_head = $x_start + $arrow_width;
    
    my @top_point    = ($x_start, $y_start);
    my @bottom_point = ($x_start, $y_start + $height);
    my @arrow_head   = ($x_arrow_head, $y_start + $side_length);

    # my @apex = ($x_center, $y_start);
    # my @left = ($x_center - $side_length, $y_start - $side_length);
    # my @right = ($x_center + $side_length, $y_start - $side_length);
    return new Sanger::Graphics::Glyph::Poly(

        # {   'points'       => [ @left, @right, @apex ],
        {   'points'       => [ @top_point, @bottom_point, @arrow_head ],
            'colour'       => $colour,
            'absolutex'    => $absolutex,
            'absolutey'    => $absolutey,
            # 'bordercolour' => 'black',
        }
    );
}

sub unique_slice_attribute {
    my ($slice, $attribute_name) = @_;
    my $value = undef;
    for my $attrib (@{ $slice->get_all_Attributes($attribute_name) }) {
        $value = $attrib->value;
    }
    return $value;
}

1;
