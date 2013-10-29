package Bio::EnsEMBL::GlyphSet::improved_sequence;
use strict;
use vars qw(@ISA);
use Bio::EnsEMBL::Feature;
use Bio::EnsEMBL::GlyphSet_simple;
use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::GlyphSet_simple);

use Readonly;

Readonly my $ODD_SCAFFOLD  => 1;
Readonly my $EVEN_SCAFFOLD => 2;
Readonly my $NO_SCAFFOLD   => 3;

Readonly my %CONTIG_COLOURS => +{
    $ODD_SCAFFOLD  => [ 'azure1',     'azure2' ],
    $EVEN_SCAFFOLD => [ 'seashell1',  'seashell2' ],
    $NO_SCAFFOLD   => [ 'whitesmoke', 'gray90' ],
};

Readonly my $GLYPH_COLOUR => 'forestgreen';
Readonly my $TAG_COLOUR   => 'darkseagreen1';

my %colour_cache = ();

sub squish { 1; }

sub my_label {
    my $self = shift;
    if ($self->{'config'}->get('lab') eq 'off') {
        return;
    }
    return 'Improved Regions';
}

sub features {
    my ($self) = @_;

    # Improved sequences are hung on the contigs, so we need to
    # (a) Project the clone slice to individual contig slices
    # (b) Get all improved_sequence features from the contigs
    # (c) Transform the sequences back onto the clone slice
    #
    # IMPORTANT! improved_sequence features should remain hung on the
    # contigs so that if contigs are shuffled on the clone, these won't
    # necessarily be affected.
    my @features = ();
    for my $contig_projection (@{ $self->{'container'}->project('contig') }) {
        my $slice          = $contig_projection->to_Slice;
        # my $contig_feature = $self->contig_feature_glyph($contig_projection);
        # push @features, $contig_feature;

        my @improved_regions = map { $_->transfer($self->{'container'}) }
            @{ $slice->get_all_MiscFeatures('improved_sequence') };
        push @features, @improved_regions;
    }
    return \@features;
}

sub tag {
    my ($self, $f) = @_;
    # return unless $f->{'_tag_colour'};
    # my $colour = $f->{'_tag_colour'};
    return {
        'style'  => 'join',
        'tag'    => $f->{'start'} . '-' . $f->{'end'},
        'colour' => $TAG_COLOUR,
        'zindex' => $f->{'z'} || -90
    };
}

sub image_label {
    return undef;
}

sub colour {
    my ($self, $f) = @_;
    return $f->{'_glyph_colour'} || $GLYPH_COLOUR;
}

sub contig_feature_glyph {
    my $self = shift;
    my ($segment) = @_;

    my $Container = $self->{'container'};
    my $start     = $segment->from_start;
    my $end       = $segment->from_end;
    my $ctg_slice = $segment->to_Slice;
    my $ORI       = $ctg_slice->strand;
    my $feature   = new Bio::EnsEMBL::Feature(
        -START  => $start,
        -END    => $end,
        -STRAND => $ORI,
    );

    my $contig_color_index = $NO_SCAFFOLD;
    $self->{'colour_position'} ||= 0;
    for my $scaffold (@{ $ctg_slice->get_all_Attributes('contig-scaffold') })
    {
        my $scaffold_name = $scaffold->value;
        my ($scaffold_number) = ($scaffold->value =~ m/Scaffold(\d+)/);
        if (defined $scaffold_number) {
            if ($scaffold_number % 2 == 0) {
                $contig_color_index = $EVEN_SCAFFOLD;
            } else {
                $contig_color_index = $ODD_SCAFFOLD;
            }
        }
    }
    my $colour = $colour_cache{ $ctg_slice->get_seq_region_id }
        ||= $CONTIG_COLOURS{$contig_color_index}
        [ $self->{'colour_position'} ];
    $feature->{'_tag_colour'}   = $colour;
    $feature->{'_glyph_colour'} = 'transparent';
    $feature->{'z'}             = -21;

    $self->{'colour_position'} = ($self->{'colour_position'} + 1) % 2;

    return $feature;
}

1;
