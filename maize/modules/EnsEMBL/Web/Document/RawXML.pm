package EnsEMBL::Web::Document::RawXML;

use strict;
use EnsEMBL::Web::Document::Common;

our @ISA = qw(EnsEMBL::Web::Document::Common);

use Data::Dumper qw(Dumper);

sub _initialize_HTML {
  my $self = shift;

## General layout for popup pages...

  $self->add_head_elements qw(
    title      EnsEMBL::Web::Document::HTML::Null
    stylesheet EnsEMBL::Web::Document::HTML::Null
    javascript EnsEMBL::Web::Document::HTML::Null
  );

  $self->add_body_elements(
    masthead       => 'EnsEMBL::Web::Document::HTML::Null',
    release        => 'EnsEMBL::Web::Document::HTML::Null',
    helplink       => 'EnsEMBL::Web::Document::HTML::Null',
    menu           => 'EnsEMBL::Web::Document::HTML::Null',
    content        => 'EnsEMBL::Web::Document::HTML::Content',
    copyright      => 'EnsEMBL::Web::Document::HTML::Null',
  );
  #$self->call_child_functions( 'common_page_elements' );
  $self->_common_HTML;
  $self->_script_HTML;
}

1;
