=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016] EMBL-European Bioinformatics Institute

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

package EnsEMBL::Web::Component::Gene::SorghumbaseSearch;

### Link from a Sorghum bicolor BTx623 v3 gene page to the same gene on SorghumBase

use strict;

use base qw(EnsEMBL::Web::Component::Gene);

sub _init {
  my $self = shift;
  $self->cacheable(0);
  $self->ajaxable(1);
}

sub content {
  my $self           = shift;
  my $gene_stable_id = $self->object->gene->stable_id;

  my $grm_server = 'sorghumbase.org/genes';
  # no '/' before '?': sorghumbase.org/genes/?idList= returns 404
  my $html .= '<br><p><a href="https://'.$grm_server. '?idList=' . $gene_stable_id . '">' . $gene_stable_id . '</a></p>';

  return $html;
}

1;
