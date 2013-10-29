package Bio::EnsEMBL::GlyphSet::marker_overgo;
use strict;
use vars qw(@ISA);
use Bio::EnsEMBL::GlyphSet_feature;

@ISA = qw(Bio::EnsEMBL::GlyphSet_feature);

sub my_label { return "Overgo Markers"; }

sub features {
    my ($self) = @_;
    my $slice = $self->{'container'};

    return $slice->get_all_MarkerFeatures('OVERGO');

}

sub href {
    my( $self, $id ) = @_;
    return $self->ID_URL( 'OVERGO', $id );
}

sub zmenu {
    my( $self, $id ) = @_;
    return { 'caption' => "OVERGO ".$id, "$id" => $self->href( $id ) };
}
1;
