package EnsEMBL::Maize::Configuration::Translation;
use strict;
use EnsEMBL::Web::Configuration::Translation;
our @ISA = qw( EnsEMBL::Web::Configuration::Translation );

# Maize-specific TransView configuration.
sub protview{
  my $self   = shift;

  my $document = $self->{page};
  my $content  = $document->content;

  # Top (Translation) Panel
  my( $panel ) = grep /^info/, $content->panels;
  $panel = $content->panel($panel);

  $panel->replace_component
      qw( name EnsEMBL::Maize::Component::Gene::name );
  $panel->remove_component
      qw( stable_id );
  $panel->remove_component
      qw( method );
  $panel->add_component_before
      qw( interpro go EnsEMBL::Maize::Component::Transcript::go );

  # Remove DAS Panel
  my( $das_panel ) = grep /^dasinfo/, $content->panels;
  $content->remove_panel( $das_panel );
}

sub context_menu{
  # Needed to prevent duplicated menu
}
