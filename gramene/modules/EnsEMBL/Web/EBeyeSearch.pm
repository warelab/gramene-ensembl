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

package EnsEMBL::Web::EBeyeSearch;

### Gramene overrides of eg-web-search's species_path and filter_species.
### Search hits name the species by production name (e.g. avena_sativa_gca947310285v1cm),
### but some genomes use a different URL name here (Avena_sativa_aslak). The original
### treated those as foreign species and linked them to http://plants.ensembl.org/...

use strict;

sub filter_species {
  ## The species filter dropdown sends the production name, but the autocomplete box
  ## (used when >200 species match) sends ucfirst'd names like Zea_mays. Gramene's
  ## search backend matches system_name case-sensitively, so convert to production name.
  my $self    = shift;
  my $species = $self->hub->param('filter_species') or return;

  return $self->hub->species_defs->get_config($species, 'SPECIES_PRODUCTION_NAME') || lc $species;
}

sub species_path {
  my ($self, $species, $genomic_unit, $want_ensembl) = @_;
  my $species_defs = $self->hub->species_defs;

  # map production name to this site's species URL name when we host it
  my $url_name = $species_defs->prodname_to_url($species);
  $species     = $url_name if $species_defs->valid_species($url_name);

  my $path = $species_defs->species_path(ucfirst($species));

  if ($path =~ /^\/$species/i and !$species_defs->valid_species(ucfirst $species) and $genomic_unit) {
    # there was no direct mapping in current unit, use the genomic_unit to add the subdomin
    $path = sprintf 'http://%s.ensembl.org/%s', $genomic_unit, $species;
  }

  # If species is in both Ensembl and EG, then $species_defs->species_path will
  # return EG url by default - sometimes we know we want ensembl
  $path =~ s/http:\/\/[a-z]+\./http:\/\/www\./ if $want_ensembl;

  return $path;
}

1;
