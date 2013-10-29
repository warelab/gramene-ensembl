package EnsEMBL::Web::Document::HTML::Null;
use strict;

use EnsEMBL::Web::Document::HTML;
use EnsEMBL::Web::Document::HTML::MastHead;
use EnsEMBL::Web::Document::HTML::SearchBox;
use EnsEMBL::Web::Document::HTML::Content;
use EnsEMBL::Web::Document::HTML::Copyright;
use EnsEMBL::Web::Document::HTML::Menu;
use EnsEMBL::Web::Document::HTML::Release;
use EnsEMBL::Web::Document::HTML::HelpLink;
use EnsEMBL::Web::Document::HTML::Title;
use EnsEMBL::Web::Document::HTML::Stylesheet;
use EnsEMBL::Web::Document::HTML::Javascript;
#use EnsEMBL::Web::Document::HTML::RSS;
#use EnsEMBL::Web::Document::HTML::Metax;

use vars qw(@ISA);
@ISA = qw(EnsEMBL::Web::Document::HTML
          EnsEMBL::Web::Document::HTML::MastHead
          EnsEMBL::Web::Document::HTML::SearchBox
          EnsEMBL::Web::Document::HTML::Content
          EnsEMBL::Web::Document::HTML::Copyright
          EnsEMBL::Web::Document::HTML::Menu
          EnsEMBL::Web::Document::HTML::Release
          EnsEMBL::Web::Document::HTML::HelpLink
          EnsEMBL::Web::Document::HTML::Title
          EnsEMBL::Web::Document::HTML::Stylesheet
          EnsEMBL::Web::Document::HTML::Javascript
          EnsEMBL::Web::Document::HTML::RSS
          EnsEMBL::Web::Document::HTML::Metax
          );

sub render {
  # Do nothing
}

1;

