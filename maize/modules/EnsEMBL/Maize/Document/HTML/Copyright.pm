package EnsEMBL::Maize::Document::HTML::Copyright;
use strict;
use CGI qw(escapeHTML);
use EnsEMBL::Web::Document::HTML;
use Maize::Page;
use SiteDefs;

our @ISA = qw(EnsEMBL::Web::Document::HTML);

sub render {
  my $self = shift;

  # We have to atach some variables to the apache request as these are
  # required by Maize::Page. The variables themselves are set in the
  # SiteDefs.pm plugin.
  my $page = Maize::Page->new() || 
      die "Could not create a Maize::Page";
  my $footer = $page->end_body;

  $self->print($footer);
}

1;

