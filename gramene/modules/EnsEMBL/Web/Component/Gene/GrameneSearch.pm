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

package EnsEMBL::Web::Component::Gene::GrameneSearch;

use strict;

use EnsEMBL::Web::Document::Image::R2R;

use base qw(EnsEMBL::Web::Component::Gene);

sub _init {
  my $self = shift;
  $self->cacheable(0);
  $self->ajaxable(1);
}

sub content {
  my $self          = shift;
  my $hub           = $self->hub;
  my $object        = $self->object;
  my $gene          = $object->gene;
  my $species_defs  = $hub->species_defs;
  my $table         = $self->new_twocol;
  my $site_type     = $species_defs->ENSEMBL_SITETYPE;
  my @CCDS          = @{$object->Obj->get_all_DBLinks('CCDS')};
  my @Uniprot       = @{$object->Obj->get_all_DBLinks('Uniprot/SWISSPROT')};
  my $db            = $object->get_db;
  my $alt_genes     = $self->get_matches('alternative_genes', 'Alternative Genes', 'ALT_GENE', 'show_version'); #gets all xrefs, sorts them and stores them on the object. Returns HTML only for ALT_GENES
  my @RefSeqMatches = @{$gene->get_all_Attributes('refseq_compare')};
  my $display_xref  = $gene->display_xref;
  my ($link_url)    = $display_xref ? $self->get_gene_display_link($gene, $display_xref) : ();

  my $gene_stable_id = $gene->stable_id;
  #my $sbase_link = '<a href="https://oryza.gramene.org/?idList=' . $gene_stable_id . '"><strong>OryzaGramene Genes search</strong></a>';

#  $table->add_row('Gramene Search', $sbase_link);

  #my $html = '<a href="https://oryza.gramene.org/?idList=' . $gene_stable_id . '"><img alt="GrameneSearchExample" src="/i/48/grmsearch.png" width="150" hight="150"></a>'; 

  my $grm_server = $species_defs->GRM_SERVERNAME;
  #my $html = '<img alt="GrameneSearchExample" src="/i/48/grmsearch.png" width="150" hight="150">';
  #https://gramene.org/?fq_field=text&fq_value=NAC13&category=text&name=NAC13
  my $html .= '<br><p><a href="https://'.$grm_server. '/?fq_field=id&fq_value=' . $gene_stable_id . '&category=Gene&name=' . $gene_stable_id . '">' .$gene_stable_id . '</a></p>';


  return $html;

}


1;
