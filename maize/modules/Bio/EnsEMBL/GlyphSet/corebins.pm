package Bio::EnsEMBL::GlyphSet::corebins;
use strict;
use vars qw(@ISA);
#use Bio::EnsEMBL::GlyphSet_feature;
use Bio::EnsEMBL::GlyphSet_simple;

use Data::Dumper qw(Dumper);    # For debug


#@ISA = qw(Bio::EnsEMBL::GlyphSet_feature);
@ISA = qw(Bio::EnsEMBL::GlyphSet_simple);

sub my_label { return "Virtual Core Bins"; }

sub features {
    my ($self) = @_;
    my $slice = $self->{'container'};

    return $slice->get_all_MiscFeatures('core_bins');

}

sub image_label {
    my ($self, $f) = @_;
    return (@{ [ $f->get_scalar_attribute('name') ] }, 'overlaid');
}

#sub zmenu {
#    my ($self, $f) = @_;
#    my $name = $f->get_scalar_attribute('name');

#    my $zmenu = {
#        'caption' => "FPC contig: $name",
#        "01:bp: @{[$f->seq_region_start]}-@{[$f->seq_region_end]}" => '',
#        "02:length: @{[$f->length]} bps"                           => '',
#        "03:Centre on FPC ctg" => $self->href($f),
#    };

#    return $zmenu;
#}

sub href {
    my ($self, $f) = @_;
    return
        "/@{[$self->{container}{_config_file_name_}]}/$ENV{'ENSEMBL_SCRIPT'}?map\
frag=@{[$f->get_scalar_attribute('name')]}";
}


#sub image_label {
#    (my $self, my $id) = @_;
#    warn Dumper ($self);
#    return ("@{[$f->get_scalar_attribute('name')]}", 'overlaid');


#    return ($f -> id, 'overlaid');
#}

#sub href {
#    my( $self, $id ) = @_;
#    return $self->ID_URL( 'Core Bin', $id );
#}

#sub zmenu {
#    my( $self, $id ) = @_;
#    return { 'caption' => "Core Bin ".$id, "$id" => $self->href( $id ) };
#}
1;
