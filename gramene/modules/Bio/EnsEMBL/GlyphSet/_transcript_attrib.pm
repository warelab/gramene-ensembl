package Bio::EnsEMBL::GlyphSet::_transcript_attrib;

use strict;
use base qw(Bio::EnsEMBL::GlyphSet::_transcript);

sub features {
    my ($self) = @_;
    my $slice = $self->{'container'};

    my $db_alias    = $self->my_config('db');
    my $attrib      = $self->my_config('attrib_code');
    # my $value       = $self->my_config('attrib_value');
    my $transcripts = $slice->get_all_Genes_by_attribute($attrib);
    $self->timer_push('Fetched transcripts', undef, 'fetch');
    $self->{geneset} = $attrib;
    return $transcripts;
}

sub gene_key {
    my ($self, $gene) = @_;
    my $value = $self->attribute_value($gene, $self->{geneset});
    if ($value ne "1") {
        return lc "$self->{geneset}-$value";
    }
    return lc $self->attribute_value($gene, 'peptide-class');
}

sub transcript_key {
    my ($self, $transcript, $gene) = @_;
    my $value = $self->attribute_value($gene, $self->{geneset});
    if ($value ne "1") {
        return lc "$self->{geneset}-$value";
    }
    return lc $self->attribute_value($transcript, 'peptide-class');
}

sub attribute_value {
    my $self = shift;
    my ($object, $attrib_code) = @_;
    my $value = "";
    eval {
        $value = $object->get_all_Attributes($attrib_code)->[0]->value();
    };
    return $value;
}

sub text_label {
    my $self = shift;
    my ($gene, $transcript) = @_;
    my @result;
    my @labels
        = split("\n", $self->SUPER::text_label($transcript, $gene));
    my $sel_t     = $self->core('t');
    my $sel_g     = $self->my_config('g') || $self->core('g');
    my $decorator = sub {
        my ($sel, $obj) = @_;
        my $label = shift @labels;
        push @result, $label;
        if ($sel eq $obj->stable_id && $label ne $sel) {
            push @result, "($sel)";
        }
    };
    if (scalar @labels == 2) {
        $decorator->($sel_t, $transcript);        
    }
    $decorator->($sel_g, $gene);
    return join("\n", @result);
}

1;

package Bio::EnsEMBL::Slice;

use Bio::PrimarySeqI;

use base qw(Bio::PrimarySeqI);

sub get_all_Genes_by_attribute {
    my $self = shift;
    my ($attribute) = @_;

    if (!$self->adaptor()) {
        warning('Cannot get Genes without attached adaptor');
        return [];
    }

    my $adaptor = $self->adaptor->db->get_GeneAdaptor();

    my $genes = $adaptor->fetch_all_by_Slice($self);

    return [ grep { scalar(@{ $_->get_all_Attributes($attribute) }) > 0 }
            @$genes ];
}

1;
