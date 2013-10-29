package EnsEMBL::Web::Document::Element::ToolLinks;

### Generates links to site tools - BLAST, help, login, etc (currently in masthead)

use strict;
use EnsEMBL::Web::Document::Element;
use CGI qw(escape);

our @ISA = qw(EnsEMBL::Web::Document::Element);

sub content   {
  my $self    = shift;
  my $dir = '/'.$ENV{'ENSEMBL_SPECIES'};
  $dir = '' if $dir !~ /_/;
  my $url     = CGI::escape($ENV{'REQUEST_URI'});
  my $html = '<div class="print_hide">';

  my $blast_dir='Multi';
  
  my $sp_param = $ENV{'ENSEMBL_SPECIES'} =~ /_/ ? "?species=$ENV{'ENSEMBL_SPECIES'}" : '';
  $html .= qq(<a href="/$blast_dir/blastview$sp_param">BLAST</a> &nbsp;|&nbsp;) if $self->blast;
  $html .= qq(<a href="/biomart/martview">BioMart</a> &nbsp;|&nbsp;)   if $self->biomart;
  $html .= qq(<a href="/tools.html">Tools</a> &nbsp;|&nbsp;) ;
  $html .= qq(<a href="/info/docs/index.html" id="help">Documentation</a> &nbsp;|&nbsp;);
  $html .= qq(<a href="/db/help" id="help">Help</a> &nbsp;|&nbsp;);
  $html .= qq(<a href="/db/feedback/send_feedback">Feedback</a>);
  $html .= '</div>';

  return ($html);
}

1;

