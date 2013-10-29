package Bio::EnsEMBL::GlyphSet::bac_map;

use strict;
use vars qw(@ISA);
use Bio::EnsEMBL::GlyphSet_simple;
@ISA = qw(Bio::EnsEMBL::GlyphSet_simple);

use Readonly;

use constant ATTRIBUTE_IS_EXTERNAL      => 'external';
use constant ATTRIBUTE_ANNOTATION_LEVEL => 'annotation';
use constant ATTRIBUTE_STATE            => 'state';
use constant ATTRIBUTE_STATUS           => 'seqstatus';
use constant ATTRIBUTE_HIGH_CONFIDENCE  => 'hc_marker';

Readonly my $TRUE  => 1;
Readonly my $FALSE => 0;

Readonly my %CLONE_CLASS_FOR => (
    'I'   => 'Genbank Record',
    'II'  => 'Primary Annotation',
    'III' => 'Secondary Annotation',
);

Readonly my @CLONE_CLASSES =>
    ('FULLTOP', 'PREFIN', 'ACTIVEFIN', 'IMPROVED', 'External');

sub my_label { return "FPC Assembly"; }

## Retrieve all BAC map clones - these are the clones in the
## subset "bac_map" - if we are looking at a long segment then we only
## retrieve accessioned clones ("acc_bac_map")

sub features {
    my ($self) = @_;
    my @sorted
        = sort { $a->seq_region_start <=> $b->seq_region_start }
        @{ $self->{'container'}->get_all_MiscFeatures($self->get_set_name) };
    my @clone_legend
        = map { ($_ => $self->{'colours'}{"col_${_}"}) } @CLONE_CLASSES;
    $self->{'config'}->{'clone_legend_features'}->{'clones'} = +{
        'priority' => 1000,
        'legend'   => \@clone_legend,
    };
    return \@sorted;
}

=pod

=head2 get_set_name
    Returns the name of the feature set that contains the desired features

=cut

sub get_set_name {
    return 'bac_map';
}

## If bac map clones are very long then we draw them as "outlines" as
## we aren't convinced on their quality...

sub colour {
    my $self = shift;
    my ($f) = @_;

    my $color_key    = undef;
    my $border       = undef;
    my $label_colour = $self->{'colours'}{'label'};

    my $status   = $f->get_scalar_attribute(ATTRIBUTE_STATUS);
    my $external = $f->get_scalar_attribute(ATTRIBUTE_IS_EXTERNAL);
    my $annotation_level
        = $f->get_scalar_attribute(ATTRIBUTE_ANNOTATION_LEVEL);
    if ($external eq 'true') {
        $color_key = 'col_External';
        $label_colour = 'white';
    } else {
        $color_key = "col_${status}";
        $label_colour = 'white' if length($status) > 0;
    }
    return $self->{'colours'}{$color_key}, $label_colour, $border;
}

## Return the image label and the position of the label
## (overlaid means that it is placed in the centre of the
## feature.

sub image_label {
    my ($self, $f) = @_;
    return ("@{[$f->get_scalar_attribute('name')]}", 'overlaid');
}

## Link back to this page centred on the map fragment

sub href {
    my ($self, $f) = @_;
    return
        "/@{[$self->{container}{_config_file_name_}]}/$ENV{'ENSEMBL_SCRIPT'}?mapfrag=@{[$f->get_scalar_attribute('name')]}";
}

sub tag {
    my ($self, $f) = @_;
    my @result = ();
    my $bef    = $f->get_scalar_attribute('BACend_flag');

    my @bacends = @{ $f->get_all_attribute_values('bacend') };
    my @ssrs    = @{ $f->get_all_attribute_values('bacend_ssr') };

    if ($self->is_hybridized_marker($f)) {
        push @result,
            +{
            'style'  => 'rect',
            'colour' => 'yellow',
            'start'  => $f->start,
            'end'    => $f->end,
            };
    }

    if    (@bacends == 2) { $bef = 3 }
    elsif (@bacends == 1) { $bef = 1 }

    (my $state = $f->get_scalar_attribute(ATTRIBUTE_STATE)) =~ s/^\d\d://;
    my ($s, $e) = $self->sr2slice(
        $f->get_scalar_attribute('inner_start'),
        $f->get_scalar_attribute('inner_end')
    );

    push @result,
        +{
        'style'  => 'right-end',
        'colour' => (scalar(@ssrs) ? 'red' : $self->{'colours'}{"bacend"}),
        }
        if ($bef == 2 || $bef == 3);
    push @result,
        +{
        'style'  => 'left-end',
        'colour' => (scalar(@ssrs) ? 'red' : $self->{'colours'}{"bacend"}),
        }
        if ($bef == 1 || $bef == 3);

    my $fp_size = $f->get_scalar_attribute('fp_size');
    if ($fp_size && $fp_size > 0) {
        my $start = int(($f->start + $f->end - $fp_size) / 2);
        my $end   = $start + $fp_size - 1;
    }
    return @result;
}

sub is_hybridized_marker {
    my $self = shift;
    my ($feature) = @_;

    my $return_value = $FALSE;

    my $feature_markers = $feature->get_all_attribute_values('clone_marker');

    my %parameters         = Apache2::RequestUtil->request->args();
    my $highlighted_marker = $parameters{'hl_marker'};

    if (   defined($highlighted_marker)
        && $highlighted_marker ne ''
        && scalar(grep {m/^$highlighted_marker$/} @{$feature_markers}) > 0)
    {
        $return_value = $TRUE;
    }
    return $return_value;
}

## Create the zmenu...
## Include each accession id separately

sub zmenu {
    my ($self, $f) = @_;
    return
        if $self->{'container'}->length() > (
                $self->{'config'}->get($self->check(), 'threshold_navigation')
                    || 2e7
        ) * 1000;

    my $name = $f->get_scalar_attribute('name');
    my $zmenu = { "caption" => "Clone: $name" };

    my $i = 0;

    $zmenu->{ sprintf("%2.2d:Centre on clone", ++$i) } = $self->href($f);

    my @accessions       = @{ $f->get_all_attribute_values('embl_acc') };
    my $annotation_level = $f->get_scalar_attribute('annotation');

    if ($annotation_level ne 'I') {
        my $clone_species = 'Zea_mays2';
        if ($f->get_scalar_attribute('external') eq 'true') {
            $clone_species = 'Zea_mays_external';
        }
        for (@accessions) {
            $zmenu->{ sprintf("%2.2d:Jump to BAC View", ++$i) }
                = "/$clone_species/contigview?contig=$_";    #GRAMENE
        }
    } else {
        $zmenu->{ sprintf("%2.2d:BAC View not available", ++$i) } = '';
    }

    for (@accessions) {
        $zmenu->{ sprintf("%2.2d:Accession: %s", ++$i, $_) }
            = "r?d=ENTREZ_NUCLEOTIDE&ID=$_";
    }
    for (@{ $f->get_all_attribute_values('bacend') }) {      # GRAMENE
        $zmenu->{ sprintf("%2.2d:BACend: $_", ++$i) }
            = "r?d=ENTREZ_NUCLEOTIDE&ID=$_";
    }
    for (@{ $f->get_all_attribute_values('bacend_ssr') }) {    # GRAMENE
        $zmenu->{ sprintf("%2.2d:BACend SSR: $_", ++$i) } = '';
    }

    for (@{ $f->get_all_attribute_values('clone-overgo') }) {
        $zmenu->{ sprintf("%2.2d:Overgo: $_", ++$i) } = '';
    }

    # Extra zmenu from EXTURL hack
    if (my $exturls = $self->species_defs->ENSEMBL_EXTERNAL_URLS) {
        if (my $label = $exturls->{EXT_CLONE_URL_LABEL}) {
            $zmenu->{ sprintf("%2.2d:$label: $name", ++$i, $name) }
                = "r?d=EXT_CLONE_URL&ID=$name";
        }
    }

    # Extra zmenus can be configured
    my %extz = (
        %{ $self->{'config'}->get("bac_map", 'ZMENU') || {} },
        %{ $self->my_config('ZMENU') || {} }
    );
    for (keys %extz) {
        $zmenu->{ sprintf("%2.2d: %s", ++$i, $_) } = $extz{$_};
    }

    my $start  = $f->seq_region_start;
    my $end    = $f->seq_region_end;
    my $length = $f->length;
    $zmenu->{ sprintf("%2.2d:bp: %d-%d", ++$i, $start, $end) } = '';
    $zmenu->{ sprintf("%2.2d:length: %d", ++$i, $length) } = '';

    (my $state = $f->get_scalar_attribute('state')) =~ s/^\d\d://;
    my $bac_info = ('Interpolated', 'Start located', 'End located',
        'Both ends located')[ $f->get_scalar_attribute('BACend_flag') ];

    if (my $org = $f->get_scalar_attribute('organisation')) {
        $zmenu->{ sprintf("%2.2d::Organisation: %s", ++$i, $org) } = '';
    }
    if ($state) {
        $zmenu->{ sprintf("%2.2d:State: %s", ++$i, $state) } = '';
    }
    if (my $len = $f->get_scalar_attribute('seq_len')) {
        $zmenu->{ sprintf("%2.2d:Seq length: %s", ++$i, $len) } = '';
    }
    if (my $status = $f->get_scalar_attribute(ATTRIBUTE_STATUS)) {
        $zmenu->{ sprintf("%2.2d:Sequence Status: %s", ++$i, $status) } = '';
    }
    if ($annotation_level) {
        $zmenu->{
            sprintf("%2.2d:Clone Class: %s",
                ++$i, $CLONE_CLASS_FOR{$annotation_level})
            }
            = '';
    }

    if (my $ctg = $f->get_scalar_attribute('superctg')) {
        $zmenu->{ sprintf("%2.2d:FPContig: %s", ++$i, $ctg) } = '';
    }

    return $zmenu;
}

sub legend {
    my ($self, $colours) = @_;
    return (
        'clones', 1011,
        [   'External' => 'peru',
            'Level I'  => 'seagreen2',
        ]
    );
}

1;
