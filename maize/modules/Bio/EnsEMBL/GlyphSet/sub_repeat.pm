package Bio::EnsEMBL::GlyphSet::sub_repeat;
use strict;
use vars qw(@ISA);
use Bio::EnsEMBL::GlyphSet::repeat_lite;
@ISA = qw( Bio::EnsEMBL::GlyphSet::repeat_lite );

sub my_label {
    my $self = shift;
    return $self->shortened_name
}

sub check { return 'sub_repeat'; }

sub repeat_type {
    my $self = shift;
    return $self->{'extras'}->{'name'};
}

sub managed_name {
    my $self = shift;
    my $managed = $self->shortened_name;
    my %roman = ( 'I' => 1, 'II' => 2, 'III' => 3);
    $managed =~ s/\b(I+)\b/$roman{$1}/g;
    $managed =~ s/\W+/_/g;
    $managed = 'Unknown' if ($managed !~ m/\w/);
    return lc "managed_repeat_$managed";
}

sub shortened_name {
    my $self = shift;
    my $name = $self->repeat_type;
    $name =~ s!/REcat!!s;
    $name =~ s!(Class \S+).*$!$1!;
    return $name;
}

1;
