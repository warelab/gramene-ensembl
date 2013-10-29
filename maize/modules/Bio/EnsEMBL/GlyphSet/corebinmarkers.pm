package Bio::EnsEMBL::GlyphSet::corebinmarkers;
use strict;
use vars qw(@ISA);
use Bio::EnsEMBL::GlyphSet_simple;

@ISA = qw(Bio::EnsEMBL::GlyphSet_simple);

#sub my_label { return "Overgo Markers"; }

sub features {
    
    my ($self) = @_;
    my $slice = $self->{'container'};
    
    warn $slice;
    warn join ", " , @{$slice->get_all_MarkerFeatures('core_bin_marker')};
    return $slice->get_all_MarkerFeatures('core_bin_marker');
    

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
