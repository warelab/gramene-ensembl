package EnsEMBL::Maize::Component::Transcript;
use strict;
use EnsEMBL::Web::Component::Transcript;
our @ISA = qw( EnsEMBL::Web::Component::Transcript);

# GRAMENE - We change the description of how the terms were mapped,
# and link to Maize ontology rather than Ensembl GoView
sub go{
  my( $panel, $object ) = @_;
  my $label = 'GO';

  # TODO: If object is a Translation, get the Transcript object
  if( $object->__objecttype eq 'Translation'){
    $object = $object->get_transcript_object;
  }

  unless ($object->__data->{'go_links'}){
    # Nothing cached. Generate similarity links.
    my @similarity_links = @{$object->get_similarity_hash($object->Obj)};
    return unless (@similarity_links);
    EnsEMBL::Web::Component::Transcript::_sort_similarity_links
        ($object, @similarity_links);
  }
  my $go_hash = $object->get_go_list() || return;

  my $html = ('<dl><dt><strong>The following GO terms have been mapped to '. 
              'this entry via '.
              '<a href="http://www.ebi.ac.uk/interpro/README1.html">'.
              '<strong>InterProScan</strong></a>:</strong></dt>' );

  foreach my $go (sort keys %{$go_hash}){
    my @go_data = @{$go_hash->{$go}||[]};
    my( $evidence, $desc ) = @go_data;
    my $link = $desc;
    $link =~ s/ /\+/g;
    my $goidurl  = $object->get_ExtURL_link($go,'GO',$go);
    my $queryurl = $object->get_ExtURL_link($desc,'GOTERMNAME',$link);
    my $evidurl  = $object->get_ExtURL_link
        ($evidence,'GOEVIDENCECODE',$evidence);
    my $tmpl = '<dd>%s [%s] <code>%s</code></dd>'."\n";
    $html .= sprintf
        ( $tmpl, $goidurl, $queryurl, $evidurl );
  }
  $html .= qq(</dl>);
  $panel->add_row( $label, $html );
}

