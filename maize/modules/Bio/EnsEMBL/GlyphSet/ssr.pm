package Bio::EnsEMBL::GlyphSet::ssr;
use strict;
use vars qw(@ISA);
use Bio::EnsEMBL::GlyphSet_simple;
@ISA = qw(Bio::EnsEMBL::GlyphSet_simple);

sub my_label { return "SSR Markers"; }

sub features {
    my ($self) = @_;
    warn("Looking for SSR features");
    my $container_length = $self->{'container'}->length();
    my $max_full_length = $self->{'config'}->get("bac_map", 'full_threshold')
        || 2e4;

    my @ssrs = $self->{'container'}->get_all_MiscFeatures('bac_map');
    warn("*** Found ", scalar @ssrs, " SSRs");
    return \@ssrs;
}

## If bac map clones are very long then we draw them as "outlines" as
## we aren't convinced on their quality...

sub colour {
    my ($self, $f) = @_;
    return $self->{'colours'}{'ssr'}, $self->{'colours'}{'label'}, undef;
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
    my @ssrs   = @{ $f->get_all_attribute_values('bacend_ssr') };

    for my $ssr (@ssrs) {
        push @result, { 'style'  => 'right-end', 'colour' => $self->{'colours'}{'ssr'} };
    }

    return @result;
}
## Create the zmenu...
## Include each accession id separately

sub zmenu {
    my ($self, $f) = @_;
    return
        if $self->{'container'}->length()
        > ($self->{'config'}->get($self->check(), 'threshold_navigation')
            || 2e7) * 1000;

    my $name  = $f->get_scalar_attribute('name');
    my $zmenu = { "caption" => "Clone: $name" };

    my $i = 0;

    $zmenu->{ sprintf("%2.2d:Centre on clone", ++$i) } = $self->href($f);

    my @accessions = @{ $f->get_all_attribute_values('embl_acc') };

    foreach (@accessions) {
        $zmenu->{ sprintf("%2.2d:Jump to BAC Browser", ++$i) }
            = "/Zea_mays2/contigview?contig=$_";    #GRAMENE
    }

    #$zmenu->{sprintf("%2.2d:Jump to Gramene Marker",++$i)} =
    #    "r?d=MARKER&ID=$name"; # GRAMENE disable for build 17a
    if (my $ctg = $f->get_scalar_attribute('superctg')) {
        $zmenu->{ sprintf("%2.2d:Jump to CMap", ++$i) }
            = "r?d=CMAP_FPC_VIEWER&ID=$ctg&ID2=$name";
    }

    foreach (@accessions) {
        $zmenu->{ sprintf("%2.2d:Accession: %s", ++$i, $_) }
            = "r?d=ENTREZ_NUCLEOTIDE&ID=$_";
    }
    foreach (@{ $f->get_all_attribute_values('bacend') }) {    # GRAMENE
        $zmenu->{ sprintf("%2.2d:BACend: $_", ++$i) }
            = "r?d=ENTREZ_NUCLEOTIDE&ID=$_";
    }
    foreach (@{ $f->get_all_attribute_values('bacend_ssr') }) {    # GRAMENE
        $zmenu->{ sprintf("%2.2d:BACend SSR: $_", ++$i) } = '';
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
    foreach (keys %extz) {
        $zmenu->{ sprintf("%2.2d: %s", ++$i, $_) } = $extz{$_};
    }

    my $start  = $f->seq_region_start;
    my $end    = $f->seq_region_end;
    my $length = $f->length;
    $zmenu->{ sprintf("%2.2d:bp: %d-%d", ++$i, $start, $end) } = '';
    $zmenu->{ sprintf("%2.2d:length: %d", ++$i, $length) } = '';

    (my $state = $f->get_scalar_attribute('state')) =~ s/^\d\d://;
    my $bac_info
        = ('Interpolated', 'Start located', 'End located', 'Both ends located')
        [ $f->get_scalar_attribute('BACend_flag') ];

    return $zmenu;
}

1;
