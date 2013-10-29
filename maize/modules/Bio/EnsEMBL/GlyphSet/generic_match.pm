package Bio::EnsEMBL::GlyphSet::generic_match;
use strict;
use Bio::EnsEMBL::GlyphSet_feature;
@Bio::EnsEMBL::GlyphSet::generic_match::ISA
    = qw(Bio::EnsEMBL::GlyphSet_feature);

use CGI;
use EnsEMBL::Maize::Util::GlyphHelper;

sub my_label {
    my $self = shift;
    return $self->my_config('TEXT_LABEL') || 'Missing label';
}

sub colour {
    my ($self, $id, $feature) = @_;
    my $colour = undef;
    if (defined $id && defined $self->object_type($id)) {
        $colour ||= $self->{'colours'}->{$self->object_type($id)};
    }
    if (defined $feature &&
        defined $feature->analysis() &&
        defined $feature->analysis()->logic_name()) {
        $colour ||= $self->{'colours'}->{$feature->analysis()->logic_name()};
    }
    $colour ||= $self->my_config('col');
    $colour ||= 'black';
    return $colour;
}

sub features {
    my ($self) = @_;

    my $method   = $self->my_config('CALL')     || 'get_all_DnaAlignFeatures';
    my $database = $self->my_config('DATABASE') || undef;
    my $threshold = (
        defined($self->my_config('THRESHOLD'))
        ? $self->my_config('THRESHOLD')
        : 80
    );
    my @logic_names;
    if ($self->my_config('FEATURES') eq 'UNDEF') {
        @logic_names = (undef());
    } else {
        @logic_names
            = split(/\s+/, ($self->my_config('FEATURES') || $self->check()));
    }
    my @feats;
    foreach my $nm (@logic_names) {
        push(@feats,
            @{ $self->{'container'}->$method($nm, $threshold, $database) });
    }
    my $colours = $self->my_config('colours');
    for my $key (sort keys %$colours) {
        EnsEMBL::Maize::Util::GlyphHelper->add_to_legend($self, $key,
            $colours->{$key});
    }

    return [@feats];
}

sub object_type {
    my ($self, $id) = @_;
    my $F = $self->my_config('SUBTYPE');
    return $self->{'type_cache'}{$id} ||= ref($F) eq 'CODE' ? $F->($id) : $F;
}

sub SUB_ID {
    my ($self, $id) = @_;
    my $T = $self->my_config('ID');
    if (ref($T) eq 'HASH') {
        $T = $T->{ $self->object_type($id) } || $T->{'default'};
    }
    return ($T && ref($T) eq 'CODE') ? &$T($id) : $id;
}

sub SUB_LABEL {
    my ($self, $id) = @_;
    my $T = $self->my_config('LABEL');
    if (ref($T) eq 'HASH') {
        $T = $T->{ $self->object_type($id) } || $T->{'default'};
    }
    return ($T && ref($T) eq 'CODE') ? &$T($id) : $id;

}

sub SUB_HREF { return href(@_); }

sub href {
    my ($self, $id) = @_;
    my $T = $self->my_config('URL_KEY');
    if (ref($T) eq 'HASH') {
        $T = $T->{ $self->object_type($id) };
    }
    return $self->ID_URL($T || 'SRS_PROTEIN', $self->SUB_ID($id));
}

sub SUB_HIGHLIGHT {
    my $self          = shift;
    my ($id)          = @_;
    my $query         = new CGI;
    my $highlight_url = $query->self_url();

    # Remove existing parameters, ensuring we are highlighting clones for
    # only one marker
    $highlight_url =~ s/hl_marker=([^;&]+)//g;
    $highlight_url =~ s/[;&]$//;

    my $extra_query = join('=', 'hl_marker', $id);
    return join(';', $highlight_url, $extra_query);
}

sub SUB_NUM_CLONES {
    my $self = shift;
    my ($id) = @_;

    return 'N/A';
}

sub zmenu {
    my ($self, $id, $feature_info, $slice_info) = @_;

    my $T = $self->my_config('ZMENU');

    if (ref($T) eq 'HASH') {
        $T = $T->{ $self->object_type($id) } || $T->{'default'};
    }
    $id =~ s/'/\'/g;    #'
    my @T = @{ $T || [] };
    my @zmenus = ('caption');
    for my $t (@T) {
        if ($t =~ m/###(\w+)###/) {
            if ($t =~ m/FEATUREVIEW/) {

                # Nasty hack. FeatureView must use $id! Always!
                $t =~ s/###(\w+)###/$self->ID_URL( $1, $id )/eg;
            } elsif ($self->can("SUB_$1")) {
                my $m = "SUB_$1";
                $t =~ s/###(\w+)###/$self->$m($id)/eg;
            } else {
                $t =~ s/###(\w+)###/$self->ID_URL($1, $self->SUB_ID($id))/eg;
            }
        }
        push @zmenus, $t;
    }
    my $zmenu = {@zmenus};

    # Gramene extra
    unless ($self->my_config('NO_ALIGNMENT')) {
        my $logic_name = ($self->my_config('FEATURES')
                || $self->my_config('LOGIC_NAME'));
        my $slice_name  = $slice_info->[0];
        my $slice_start = $slice_info->[1] - 1;
        my $chr_start   = $slice_start + ($feature_info->[0][0] || 1);
        my $chr_end     = $slice_start + ($feature_info->[-1][1] || 1);

        my $align_link
            = (   "/Multi/show_alignment.pl"
                . "?chr=$slice_name&chr_start=$chr_start&chr_end=$chr_end"
                . "&feature_id=$id"
                . "&blastdb=$logic_name.fa"
                . "&species=$ENV{ENSEMBL_SPECIES}");
    }
    return $zmenu;
}

1;
