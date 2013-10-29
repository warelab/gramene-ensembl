package Bio::EnsEMBL::GlyphSet::marker;
use strict;
use vars qw(@ISA);
use Bio::EnsEMBL::GlyphSet_simple;
use Sanger::Graphics::Glyph::Rect;
use Sanger::Graphics::Glyph::Text;
use Sanger::Graphics::Bump;

@ISA = qw(Bio::EnsEMBL::GlyphSet_simple);

my $MAP_WEIGHT = 2;
my $PRIORITY   = 50;

sub _init {
    my $self = shift;

    my $slice         = $self->{'container'};
    my $Config        = $self->{'config'};
    my $pix_per_bp    = $Config->transform->{'scalex'};
    my $bitmap_length = int($slice->length() * $pix_per_bp);
    my @bitmap;

    return unless $self->strand() == -1;
    $self->{'colours'} = $Config->get('marker', 'colours');
    my $fontname   = "Tiny";
    my $row_height = 8;
    my ($w, $h) = $Config->texthelper->px2bp($fontname);
    $w = $Config->texthelper->width($fontname);

    my $labels = $Config->get('marker', 'labels') eq 'on';

    my $priority = (
          $self->{container}{_config_file_name_} eq 'Homo_sapiens'
        ? $PRIORITY
        : undef
    );

    # GRAMENE - Hard-code logic names for quick v28-17 release fix!
    for my $f (
        @{ $slice->get_all_MarkerFeatures('Maize_marker') },
        @{ $slice->get_all_MarkerFeatures('OVERGO') },
        @{ $slice->get_all_MarkerFeatures('Locus') },
        @{ $slice->get_all_MarkerFeatures('REP') },
        @{ $slice->get_all_MarkerFeatures('STS') },
        @{ $slice->get_all_MarkerFeatures('RFLP') },
        @{ $slice->get_all_MarkerFeatures('Clone') },
        )
    {

        warn("looking at feature @{[$f->display_id]}\n");
        my $ms           = $f->marker->display_MarkerSynonym;
        my $fid          = $ms->name;
        my $type         = $self->type($f);
        my $bp_textwidth = $w * length("$fid ");
        my ($feature_colour, $label_colour, $part_to_colour)
            = $self->colour($f);
        my $href
            = "/@{[$self->{container}{_config_file_name_}]}/markerview?marker=$fid";

        my $zmenu = {
            'caption'           => "Marker: $fid",
            "01:MarkerView"     => $href,
            "02:Gramene Marker" => "r?d=MARKER&ID=$fid"
        };

        # Extra zmenu from EXTURL hack
        if (my $exturls = $self->species_defs->ENSEMBL_EXTERNAL_URLS) {
            if (my $label = $exturls->{EXT_MARKER_URL_LABEL}) {
                $zmenu->{"10:$label"} = "r?d=EXT_MARKER_URL&ID=$fid";
            }
        }

        $self->push(
            new Sanger::Graphics::Glyph::Rect(
                {   'x'         => $f->start - 1,
                    'y'         => 0,
                    'height'    => $row_height,
                    'width'     => ($f->end - $f->start + 1),
                    'colour'    => $feature_colour,
                    'absolutey' => 1,
                    'href'      => $href,

                    # GRAMENE
                    #'zmenu' => { 'caption' => ( $ms->source eq 'unists' ?
                    #                            "uniSTS:$fid" : "$type: $fid"),
                    #
                    # 'Marker info' => $href }
                    'zmenu' => $zmenu,
                }
            )
        );
        next unless $labels;
        my $glyph = new Sanger::Graphics::Glyph::Text(
            {   'x'         => $f->start() - 1,
                'y'         => $row_height + 2,
                'height'    => $Config->texthelper->height($fontname),
                'width'     => $bp_textwidth,
                'font'      => $fontname,
                'colour'    => $label_colour,
                'absolutey' => 1,
                'text'      => $fid,
                'href'      =>
                    "/@{[$self->{container}{_config_file_name_}]}/markerview?marker=$fid",
            }
        );

        my $bump_start = int($glyph->x() * $pix_per_bp);
        $bump_start = 0 if $bump_start < 0;
        my $bump_end = $bump_start + $bp_textwidth;
        next if $bump_end > $bitmap_length;
        my $row = &Sanger::Graphics::Bump::bump_row($bump_start, $bump_end,
            $bitmap_length, \@bitmap);
        $glyph->y($glyph->y() + (1.2 * $row * $h));
        $self->push($glyph);
    }
}

sub colour {
    my ($self, $f) = @_;
    my $type = $self->type($f);
    $type = lc($type);
    my $col = ();
    if ($f->analysis->logic_name eq "OVERGO") {
        $col = 'red';
    } elsif ($f->analysis->logic_name eq "Locus") {
        $col = 'blue';
    } elsif ($f->analysis->logic_name eq "STS") {
        $col = 'green';
    } elsif ($f->analysis->logic_name eq "REP") {
        $col = 'purple';
    } elsif ($f->analysis->logic_name eq "RFLP") {
        $col = 'yellow';
    } else {
        $col = $self->{'colours'}->{"$type"};
    }

# print STDERR "Color for type @{[$self->type($f)]} = $col @{[$f->analysis->logic_name]}\n";
    return $col, $col, '';
}

sub type {
    my ($self, $f) = @_;
    my $type = $f->marker->type;
    $type ||= $f->analysis->logic_name;
    return $type || '';
}

sub my_label {
    my $self = shift;
    return $self->my_config('track_label') || 'Markers';
}

1;
