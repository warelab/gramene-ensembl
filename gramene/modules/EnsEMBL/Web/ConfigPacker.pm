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

package EnsEMBL::Web::ConfigPacker;

### Gramene: drop Cactus HAL alignments (e.g. the rice 27-way) from the packed config
### unless HALXS is built. Without it, turning on one of those tracks on Location/View
### dies in GenomicAlignBlockAdaptor with "Can't locate loadable object for module HALXS".

use strict;
use warnings;
no warnings qw(uninitialized);

use previous qw(munge_databases_multi);

sub munge_databases_multi {
  my $self = shift;
  $self->PREV::munge_databases_multi(@_);

  return if eval { require Bio::EnsEMBL::Compara::HAL::HALXS::HALAdaptor; 1 };

  my $alignments = $self->db_tree->{'DATABASE_COMPARA'}{'ALIGNMENTS'} || {};
  delete $alignments->{$_} for grep { $alignments->{$_}{'type'} =~ /^CACTUS_HAL/ } keys %$alignments;
}

1;
