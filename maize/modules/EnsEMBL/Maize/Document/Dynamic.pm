package EnsEMBL::Maize::Document::Dynamic;

use strict;
use EnsEMBL::Web::Document::Common;

our @ISA = qw(EnsEMBL::Web::Document::Common);

use Data::Dumper qw(Dumper);

sub set_title {
  my $self  = shift;
  my $title = shift;
  $self->title->set( $self->species_defs->ENSEMBL_SITE_NAME.' v'.$self->species_defs->ENSEMBL_VERSION.': '.$self->species_defs->SPECIES_BIO_NAME.' '.$title );
}

sub _initialize_TextGz {
  my $self = shift; 
  $self->add_body_elements qw(
    content EnsEMBL::Web::Document::Text::Content
  );
  $self->_init();
}

sub _initialize_Text {
  my $self = shift; 
  $self->add_body_elements qw(
    content EnsEMBL::Web::Document::Text::Content
  );
  $self->_init();
}

sub _initialize_XML {
  my $self = shift; 

  $self->add_body_elements qw(
			      content     EnsEMBL::Web::Document::XML::Content
			      );

  $self->_init();
}

sub _initialize_HTML {
  my $self = shift;

## General layout for dynamic pages...

  $self->add_head_elements qw(
    title      EnsEMBL::Web::Document::HTML::Title
    stylesheet EnsEMBL::Web::Document::HTML::Stylesheet
    javascript EnsEMBL::Web::Document::HTML::Javascript
    meta       EnsEMBL::Web::Document::HTML::Meta
  );
    #iehover    EnsEMBL::Web::Document::HTML::IEHoverHack
  $self->add_body_elements qw(
    javascript_div EnsEMBL::Web::Document::HTML::JavascriptDiv
    masthead   EnsEMBL::Web::Document::HTML::MastHead
    searchbox  EnsEMBL::Web::Document::HTML::SearchBox
    release    EnsEMBL::Web::Document::HTML::Release
    helplink   EnsEMBL::Web::Document::HTML::HelpLink
    html_start EnsEMBL::Web::Document::HTML::HTML_Block
    menu       EnsEMBL::Web::Document::HTML::Menu
    content    EnsEMBL::Web::Document::HTML::Content
    copyright  EnsEMBL::Web::Document::HTML::Copyright
    html_end   EnsEMBL::Web::Document::HTML::HTML_Block
  );

  $self->call_child_functions( 'common_page_elements','static_page_elements' );

  $self->_common_HTML();
  $self->_script_HTML();
  $self->helplink->kw = $ENV{'ENSEMBL_SCRIPT'}.';se=1';
## Let us set up the search box...
  $self->searchbox->sp_common  = $self->species_defs->SPECIES_COMMON_NAME;
#  --- First the search index drop down
  if( $ENV{'ENSEMBL_SPECIES'} ne 'Multi' ) { # If we are in static content for a species
    foreach my $K ( sort @{($self->species_defs->ENSEMBL_SEARCH_IDXS)||[]} ) {
      $self->searchbox->add_index( $K );
    }
    my $T = $self->species_defs->SEARCH_LINKS || {};
    ## Now grab the default search links for the species
    my $T = $self->species_defs->SEARCH_LINKS || {};
    my $flag = 0;
    my $regexp = '^('.$ENV{'ENSEMBL_SCRIPT'}.'\d+)_URL';
    foreach my $K ( sort keys %$T ) {
      if( $K =~ /$regexp/i ) {
        $flag = 1;
        $self->searchbox->add_link( $T->{$K}, $T->{$1."_TEXT"} );
      }
    }
    unless($flag) { 
      foreach my $K ( sort keys %$T ) {
        if( $K =~ /DEFAULT(\d)_URL/ ) {
          $self->searchbox->add_link( $T->{$K}, $T->{"DEFAULT$1"."_TEXT"} );
        }
      }
    }
  } else { # If we are in general static content...
    ## Grab all the search indexes...
    foreach my $K ( $self->species_defs->all_search_indexes ) {
      $self->searchbox->add_index( $K );
    }
    ## Note we have no example links here!!
  }
#  --- and the search box links...

  $self->call_child_functions( 'extra_configuration' );
  $self->call_child_functions( 'common_menu_items', 'dynamic_menu_items' );

}

1;
