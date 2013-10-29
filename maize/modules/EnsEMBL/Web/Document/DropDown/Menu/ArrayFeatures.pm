package EnsEMBL::Web::Document::DropDown::Menu::ArrayFeatures;

use strict;
use EnsEMBL::Web::Document::DropDown::Menu;

@EnsEMBL::Web::Document::DropDown::Menu::ArrayFeatures::ISA
  =qw( EnsEMBL::Web::Document::DropDown::Menu );

sub new {
  my $class = shift;
  my $self = $class->SUPER::new(
    @_, ## This contains the menu containers as the first element
    'image_name'  => 'y-array_features',
    'image_width' => 59,
    'alt'         => 'Array-based Features'
  );
  my @menu_entries = @{$self->{'config'}->get('_settings','array_features')||[]};;
  return undef unless @menu_entries;
  foreach my $m ( @menu_entries ) {
    foreach my $c ( @{$self->{'configs'}||[]}, $self->{'config'} ) {
      if( $c->is_available_artefact($m->[0] ) ) {
        $self->add_checkbox( @$m );
        last;
      }
    }
  }
  return $self;
}

1;
