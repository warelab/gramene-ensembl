package Bio::EnsEMBL::GlyphSet::tilepath;
use strict;
use vars qw(@ISA);
use Bio::EnsEMBL::GlyphSet::bac_map;
@ISA = qw(Bio::EnsEMBL::GlyphSet::bac_map);

sub my_label { return "Sequenced BACs"; }

=pod

=head2 get_set_name
    Returns the feature set name

=cut

sub get_set_name {
    return 'acc_bac_map';
}

1;
