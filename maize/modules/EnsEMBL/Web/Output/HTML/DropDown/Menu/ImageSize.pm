package EnsEMBL::Web::Output::HTML::DropDown::Menu::ImageSize;

use strict;
use EnsEMBL::Web::Output::HTML::DropDown::Menu;
@EnsEMBL::Web::Output::HTML::DropDown::Menu::ImageSize::ISA
  =qw( EnsEMBL::Web::Output::HTML::DropDown::Menu );

sub new {
  my $class  = shift;
  my $self = $class->SUPER::new( 
    @_, ## This contains the menu containers as the first element
    'image_name'  => 'y-width',
    'image_width' => 53,
    'alt'         => 'Resize image'
  ); 
  my $LINK = sprintf qq(/%s/%s?%s), $self->{'species'}, $self->{'script'}, $self->{'LINK'};
  foreach( qw(700 900 1200 1500 2000 ) ) {
    $self->add_link( ($_==$self->{'config'}->get('_settings','width') ? "* " : '' ). "Width $_".'px', $LINK."width=$_", '' );
  }
  return $self;
}

1;
