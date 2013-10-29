package EnsEMBL::Maize::Configuration::Transcript;
use strict;
use EnsEMBL::Web::Configuration::Transcript;
our @ISA = qw( EnsEMBL::Web::Configuration::Transcript );

# Maize-specific TransView configuration.
sub transview{
  my $self   = shift;

  my $document = $self->{page};
  my $content  = $document->content;

  # Top (Transcript) Panel
  my( $panel ) = grep /^info/, $content->panels;
  $panel = $content->panel($panel);

  $panel->replace_component
      qw( name EnsEMBL::Maize::Component::Gene::name );
  $panel->remove_component
      qw( stable_id );
  $panel->remove_component
      qw( method );
  $panel->replace_component
      qw( go EnsEMBL::Maize::Component::Transcript::go );

}

sub context_menu{
  # Needed to prevent duplicated menu
}
