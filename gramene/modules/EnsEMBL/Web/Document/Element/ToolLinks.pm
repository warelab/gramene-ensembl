package EnsEMBL::Web::Document::Element::ToolLinks;

### Generates links to site tools - BLAST, help, login, etc (currently in masthead)

use strict;
use EnsEMBL::Web::Document::Element;
use CGI qw(escape);

our @ISA = qw(EnsEMBL::Web::Document::Element);

my $blastserver = 'blast.gramene.org';

sub content   {
  my $self    = shift;
  my $dir = '/'.$ENV{'ENSEMBL_SPECIES'};
  $dir = '' if $dir !~ /_/;
  my $url     = CGI::escape($ENV{'REQUEST_URI'});
  my $html = '<div class="print_hide">';

  my $blast_dir = $ENV{'ENSEMBL_SPECIES'} || 'Multi';
  
  my $sp_param = $ENV{'ENSEMBL_SPECIES'} =~ /_/ ? "?species=$ENV{'ENSEMBL_SPECIES'}" : '';
  $html .= qq(<a href="http://$blastserver/$blast_dir/blastview$sp_param">BLAST</a> &nbsp;|&nbsp;) if $self->blast;
  $html .= qq(<a href="http://ensembl.gramene.org/biomart/martview">BioMart</a> &nbsp;|&nbsp;)   if $self->biomart;
  $html .= qq(<a href="/tools.html">Tools</a> &nbsp;|&nbsp;)  ;
  $html .= qq(<a href="http://www.gramene.org/contact">Feedback</a>);
  $html .= '</div>';

  return ($html);
}

1;

