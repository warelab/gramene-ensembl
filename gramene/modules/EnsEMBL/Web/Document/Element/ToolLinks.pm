=head1 LICENSE

Copyright [2009-2014] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

package EnsEMBL::Web::Document::Element::ToolLinks;

use strict;
use warnings;
use URI::Escape;

#sub links {
#  my $self  = shift;
#  my $hub   = $self->hub;
#  my $sd    = $self->species_defs;
#  my @links;
#
#  #my $blastserver = 'blast.gramene.org';
#  #my $blast_dir =  $sd->SPECIES_SCIENTIFIC_NAME || 'Multi';
#  #used to be $ENV{'ENSEMBL_SPECIES'}, may also be $sd->SPECIES_COMMON_NAME
#
#  #my $blast_url = uri_escape("http://$blastserver/$blast_dir/blastview");
#
#  push @links, 'eqsearch',      '<a class="constant" href="/Multi/enasearch">Sequence Search</a>' if $sd->ENSEMBL_ENASEARCH_ENABLED;
#  push @links, 'blast',         '<a class="constant" href="http://blast.gramene.org/Multi/blastview">BLAST</a>' if $sd->ENSEMBL_BLAST_ENABLED;
#  push @links, 'biomart',       '<a class="constant" href="/biomart/martview">BioMart</a>';
#  push @links, 'tools',         '<a class="constant" href="/tools.html">Tools</a>';
#  push @links, 'downloads',     '<a class="constant" href="/downloads.html">Downloads</a>';
#  push @links, 'help',          '<a class="constant" href="/info/website/index.html">Help</a>';
#  push @links, 'feedback',      '<a class="constant" href="http://www.gramene.org/contact">Feedback</a>';
#  return \@links;
#}

sub links {
  my $self  = shift;
  my $hub   = $self->hub;
  my $sd    = $self->species_defs;
  my @links;

  push @links, 'eqsearch',      '<a class="constant" href="/Multi/enasearch">Sequence Search</a>' if $sd->ENSEMBL_ENASEARCH_ENABLED;

  if( $sd->ENSEMBL_BLAST_ENABLED ){
	my $blast_link = $self->hub->url({'species' => '', 'type' => 'Tools', 'action' => 'Blast'});
	$blast_link =~ s/genome_browser//i; # but genome_browser is not embedded here, it is in the base url
	push @links, 'blast', sprintf '<a class="constant" href="%s">BLAST</a>', $blast_link;
  	#push @links, 'blast', sprintf '<a class="constant" href="/%s">BLAST</a>',  if $sd->ENSEMBL_BLAST_ENABLED;
  }

  push @links, 'biomart',       '<a class="constant" href="/biomart/martview">BioMart</a>';
  push @links, 'tools',         '<a class="constant" href="/tools.html">Tools</a>';
  push @links, 'downloads',     '<a class="constant" href="/downloads.html">Downloads</a>';
  push @links, 'help',          '<a class="constant" href="/info/website/index.html">Help</a>';
#  push @links, 'docs',          '<a class="constant" href="http://www.ensemblgenomes.org/info">Documentation</a>';
  push @links, 'feedback',      '<a class="constant" href="http://dev.gramene.org/feedback">Feedback</a>'; #http://tools.gramene.org/feedback

# test upload link
# UserData/SelectFile?db=core
	
  my $upload_link = 'UserData/SelectFile?db=core';
  push @links, 'uploadData',  sprintf '<a class="constant" href="/%s">UploadData</a>', $upload_link;  
  return \@links;
}


1;

