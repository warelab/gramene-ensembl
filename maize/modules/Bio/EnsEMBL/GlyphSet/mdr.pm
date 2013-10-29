package Bio::EnsEMBL::GlyphSet::mdr;
use strict;
use vars qw(@ISA);
use Bio::EnsEMBL::GlyphSet_simple;
use Data::Dumper;
@ISA = qw(Bio::EnsEMBL::GlyphSet_simple);

sub my_label {
    my ($self) = @_;
    return $self->my_config('track_label');
}

sub features {
    my ($self) = @_;
    my $slice = $self->{'container'};

    my $logic_name = $_[0]->my_config('key');
    my @features   = @{ $slice->get_all_SimpleFeatures($logic_name) };

    return \@features;
}

sub href { my ($self, $f) = @_; return undef; }

sub zmenu {
    my ($self, $f) = @_;

    my $score = $f->can('score') ? $f->score() : '';
    my ($start, $end) = $self->slice2sr($f->start, $f->end);
    return {
        'caption'                      => $self->my_config('caption'),
        "01:Average copies:    $score" => '',
        "02:Location: $start-$end"     => '',
        "03:Length:   " . ($end - $start + 1) => ''
    };
}
1;
