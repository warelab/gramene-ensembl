package Bio::EnsEMBL::GlyphSet::supercontigs;
use strict;
use vars qw(@ISA);
use Bio::EnsEMBL::GlyphSet_simple;
@ISA = qw(Bio::EnsEMBL::GlyphSet_simple);

sub my_label { return "FPC Contigs"; }

sub features {
    my ($self) = @_;
    my $features = $self->{'container'}->get_all_MiscFeatures('superctgs');
    # local $Data::Dumper::Maxdepth = 4;
    # for my $feature (@$features) {
    #     warn Data::Dumper::Dumper($feature);
    # }
    return $features;
}

sub href {
    my ($self, $f) = @_;
    return
        "/@{[$self->{container}{_config_file_name_}]}/$ENV{'ENSEMBL_SCRIPT'}?mapfrag=@{[$f->get_scalar_attribute('name')]}";
}

sub colour {
    my ($self, $f) = @_;
    $self->{'_colour_flag'} = $self->{'_colour_flag'} == 1 ? 2 : 1;
    return $self->{'colours'}{"col$self->{'_colour_flag'}"},
        $self->{'colours'}{"lab$self->{'_colour_flag'}"};
}

sub image_label {
    my ($self, $f) = @_;
    return (@{ [ $f->get_scalar_attribute('name') ] }, 'overlaid');
}

sub zmenu {
    my ($self, $f) = @_;
    my $name = $f->get_scalar_attribute('name');

    my $zmenu = {
        'caption' => "FPC contig: $name",
        "01:bp: @{[$f->seq_region_start]}-@{[$f->seq_region_end]}" => '',
        "02:length: @{[$f->length]} bps"                           => '',
        "03:Centre on FPC ctg" => $self->href($f),
        "04:View FPC in CMap"  => "r?d=CMAP_FPC_VIEWER&ID=$name&ID2=",
    };

    # Gramene

    return $zmenu;
}

1;
