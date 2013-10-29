package EnsEMBL::Maize::Configuration::Gene;
use strict;
use EnsEMBL::Web::Configuration::Gene;
our @ISA = qw( EnsEMBL::Web::Configuration::Gene );

# Maize-specific GeneView configuration.
sub geneview{
  my $self   = shift;
  my $document = $self->{page};
  my $content  = $document->content;

  # Top (Gene) Panel
  my $panel = $content->panel('info0');
  $panel->replace_component
      qw( name EnsEMBL::Maize::Component::Gene::name );
  $panel->remove_component
      qw( stable_id );
  $panel->remove_component
      qw( diseases );

  # Remove DAS Panel
  $content->remove_panel('dasinfo0');

  # Bottom (Transcript) Panels
  PANEL: 
  foreach my $panel_name( $content->panels ){
    unless( $panel_name =~ /^trans/ ){ next PANEL }
    my $tpanel = $content->panel($panel_name);
    $tpanel->replace_component
        qw( go 
            EnsEMBL::Maize::Component::Transcript::go );
    #$tpanel->replace_component
    #    qw( similarity_matches 
    #        EnsEMBL::Maize::Component::Transcript::similarity_matches );
  }

  return 1;
}

# Swaps the name and stable_id components for Maize versions
sub _change_gene_name{
  my $self   = shift;
  my $document = $self->{page};
  my $content  = $document->content;
  # Top (Gene) Panel
  my $panel = $content->panel('info0');
  $panel->replace_component
      qw( name EnsEMBL::Maize::Component::Gene::name );
  $panel->remove_component
      qw( stable_id );

}

sub genespliceview{
  my $self   = shift;
  return $self->_change_gene_name(@_);
}

sub geneseqview{
  my $self   = shift;
  return $self->_change_gene_name(@_);
}

sub genesnpview{
  my $self   = shift;
  return $self->_change_gene_name(@_);
}

sub context_menu{
  # Needed to prevent duplicated menu
}
