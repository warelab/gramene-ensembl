package Bio::EnsEMBL::GlyphSet::repeat_lite;
use strict;
use vars qw(@ISA);
use Bio::EnsEMBL::GlyphSet_simple;
use Data::Dumper;

@ISA = qw( Bio::EnsEMBL::GlyphSet_simple );

use Readonly;

Readonly my @MIPS_COLUMNS =>
    qw(mips_recat report_class l1_class l2_class l3_class l4_class l5_class l6_class);

sub my_label {
    my $self = shift;
    return $self->my_config('track_label');
}

sub repeat_type {
    return undef;
}

sub features {
    my $self     = shift;
    my $slice    = $self->{'container'};
    my $species  = $self->my_config('species');
    my @features = sort { $a->seq_region_start <=> $b->seq_region_start }
        @{ $slice->get_all_RepeatFeatures(undef, $self->repeat_type) };

    for my $feature (@features) {
        $self->_add_mips_information($feature);
    }

    $self->{'config'}->{'repeat_legend_features'}->{'repeat_classes'} = +{
        'priority' => 1000,
        'legend'   => [
            map { $_ => $self->{colours}->{$_} }
                sort keys %{ $self->{colours} },
        ],
    };

    push @{ $self->{'config'}->{'repeat_legend_features'}->{'repeat_classes'}
            ->{'legend'} }, ('Other' => $self->my_config('col'));
    return \@features;
}

our (%mips_cache);
our $attribute_statement = undef;

sub _add_mips_information {
    my $self = shift;
    my ($feature) = @_;

    my $key = $feature->repeat_consensus->repeat_class;
    $key =~ s/\.//g;

    $feature->{_mips} = $mips_cache{$key} ||= $self->_lookup_mips($key);
}

sub _lookup_mips {
    my $self = shift;
    my ($key) = @_;

    return undef unless defined $key;

    if (!defined $attribute_statement) {
        my $species = $self->{'config'}->{'species'};
        my $adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'core');
        my $sql     = "select * from attrib_type where code = ?";
        $attribute_statement = $adaptor->dbc->prepare($sql);
    }
    $attribute_statement->bind_param(1, $key);
    $attribute_statement->execute();
    my $mips = +{};
    while (my $row = $attribute_statement->fetchrow_hashref()) {
        my @fields = split("\t", $row->{description});
        for (my $i = 0; $i < scalar @fields; $i++) {
            $mips->{ $MIPS_COLUMNS[$i] } = $fields[$i];
        }

    }

    return $mips;
}

sub colour {
    my $self = shift;
    my ($feature) = @_;

    my $colour        = undef;
    my $repeat_colour = undef;

    if (defined $feature->{_mips}) {
        my $value = $feature->{_mips}->{'l5_class'};
        $repeat_colour = $self->find_repeat_color($value);
    }
    if (defined($repeat_colour) && $repeat_colour ne '') {
        $colour = $repeat_colour;
    } else {
        $colour = $self->my_config('col');
    }
    return $colour;
}

sub find_repeat_color {
    my $self    = shift;
    my ($value) = @_;
    my $colour  = $self->{'colours'}->{ lc $value };

    return $colour;
}

sub zmenu {
    my ($self, $f) = @_;

    my ($start, $end) = $self->slice2sr($f->start(), $f->end());
    my $len = $end - $start + 1;

    my $i = 0;

    # my $index = sub { ; }
    ### Possibly should not use $f->repeat_consensus->name.... was f->{'hid'}
    my $zmenu = +{ 'caption' => $f->repeat_consensus()->name(), };

    $zmenu->{
        "@{[sprintf('%0.2d', ++$i)]}:Repeat Class: @{[$f->repeat_consensus->repeat_class]}"
        } = '';
    $zmenu->{
        "@{[sprintf('%0.2d', ++$i)]}:Repeat Type: @{[$f->repeat_consensus->repeat_type]}"
        } = '';

    if (defined $f->{_mips}) {
        $zmenu->{
            "@{[sprintf('%0.2d', ++$i)]}:Classification: $f->{_mips}->{l5_class}"
            } = '';
        $zmenu->{
            "@{[sprintf('%0.2d', ++$i)]}:RepeatMasker Class: $f->{_mips}->{l3_class}"
            } = '';
        $zmenu->{
            "@{[sprintf('%0.2d', ++$i)]}:RepeatMasker Subclass: $f->{_mips}->{l4_class}"
            } = '';
    }

    $zmenu->{"@{[sprintf('%0.2d', ++$i)]}:BP: $start-$end"} = '';
    $zmenu->{"@{[sprintf('%0.2d', ++$i)]}:Length: $len"}    = '';
    $zmenu->{"@{[sprintf('%0.2d', ++$i)]}:Lookup in MIPS"}  = $self->href($f);
    return $zmenu;
}

sub href {
    my $self = shift;
    my ($feature) = @_;
    return sprintf(
        'http://mips.gsf.de/proj/plant/webapp/recat/'
            . 'RepeatRetrival.jsp?searchcolumn=name&kw=%s&coltosort=',
        $feature->repeat_consensus->name
    );
}

1;
__END__

