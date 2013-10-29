package EnsEMBL::Maize::Component::Gene;
use strict;
use EnsEMBL::Web::Component::Gene;
our @ISA = qw( EnsEMBL::Web::Component::Gene);

sub name{
  my( $panel, $object ) = @_;
  my( $display_name, $dbname, $ext_id, $dbname_disp ) 
      = $object->display_xref();
  $display_name ||= $object->stable_id; # Revert to stable ID
  if( $ext_id ){
    $display_name 
        = $object->get_ExtURL_link( $display_name, $dbname, $ext_id );
  }
  my $display_html = sprintf('<p><strong>%s</strong>', $display_name );
  if( $dbname ){ $display_html .= sprintf( ' (%s)', $dbname ) }
  $panel->add_row( $object->type_name(), $display_html );
  return 1;
}
