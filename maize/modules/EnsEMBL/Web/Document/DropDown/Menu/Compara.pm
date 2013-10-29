package EnsEMBL::Web::Document::DropDown::Menu::Compara;

use strict;
use EnsEMBL::Web::Document::DropDown::Menu;
use Data::Dumper;
@EnsEMBL::Web::Document::DropDown::Menu::Compara::ISA
  =qw( EnsEMBL::Web::Document::DropDown::Menu );

sub new {
  my $class = shift;
  my $self = $class->SUPER::new(
    @_, ## This contains the menu containers as the first element
    'image_name'  => 'y-compara',
    'image_width' => 88,
    'alt'         => 'Compara'
  );
  my @menu_entries = @{$self->{'config'}->get('_settings','compara')||[]};

  my $num_checkboxes = 0;
  foreach ( @menu_entries ) { 
    $self->{'config'}->is_available_artefact($_->[0]) || next;
    $self->add_checkbox( @$_ );
    $num_checkboxes++;
  }
  $num_checkboxes || return;

  if( $self->{'config'}->{'compara'} ){
    my $LINK = sprintf qq(/%s/%s?%s), $ENV{'ENSEMBL_SPECIES'}, $self->{'script'}, $self->{'LINK'};
    my %species = (
      EnsWeb::species_defs->multi('BLASTZ_NET'),
      EnsWeb::species_defs->multi('BLASTZ_GROUP'),
      EnsWeb::species_defs->multi('PHUSION_BLASTN'),
      EnsWeb::species_defs->multi('BLASTZ_RECIP_NET'),
      EnsWeb::species_defs->multi('TRANSLATED_BLAT') 
    );
    foreach( keys %species ) {
      $self->add_link(
        "Add/Remove $_",
        $LINK."flip=$_",
        ''
      )
    }
  }
  return $self;
}

1;
