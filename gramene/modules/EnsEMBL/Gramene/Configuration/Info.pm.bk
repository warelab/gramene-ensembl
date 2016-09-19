package EnsEMBL::Gramene::Configuration::Info;
use strict;
use base qw( EnsEMBL::Web::Configuration );

sub populate_tree{
  # GRAMENE: Remove the What's New setting in the tree.
  my $self   = shift;
  if( my $node = $self->tree->get_node('WhatsNew') ){
    return $node->remove_node;
  }
}

1;
