package EnsEMBL::Gramene::Document::Configure;

use strict;
use base qw(EnsEMBL::Web::Root);

sub common_page_elements{
  my $self = shift;
  my $page = shift;
#  $page->remove_body_element( 'search_box' );
#  $page->replace_body_element( 'search_box', 
#                               'EnsEMBL::Gramene::Document::HTML::SearchBox' );
  return 1;
}

1;
