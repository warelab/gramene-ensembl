package Bio::EnsEMBL::GlyphSet::generic_transcript;
use strict;
use vars qw(@ISA);

use List::Util;

use Bio::EnsEMBL::GlyphSet_transcript;
use Bio::EnsEMBL::Utils::Eprof qw(eprof_start eprof_end eprof_dump);

@ISA = qw(Bio::EnsEMBL::GlyphSet_transcript);

sub analysis_logic_name {
    my $self = shift;
    return $self->my_config('LOGIC_NAME');
}

sub my_label {
    my $self        = shift;
    my $label       = $self->my_config('track_label') || 'Generic trans.';
    my $numchars    = 16;
    my $print_label = (
        length($label) > $numchars
        ? substr($label, 0, ($numchars - 2)) . '..'
        : $label
    );
    return $print_label;
}

sub colours {
    my $self = shift;
    return {
        'hi'    => $self->my_config('hi')      || 'white',
        'col'   => $self->my_config('col')     || undef,
        'super' => $self->my_config('superhi') || 'white',
        %{ $self->my_config('colours') || {} }
    };
}

sub features {
    my ($self) = @_;
    my $slice  = $self->{'container'};
    my $config = $self->{'config'};

    my @analyses = ($self->analysis_logic_name || '');

    my @biotypes = ();
    for my $biotype (
        qw/protein_coding
        protein_coding_unsupported
        transposon_pseudogene/
        )
    {
        if ($config->get('_settings', "opt_$biotype")) {
            push @biotypes, $biotype;
        }
    }
    my @genes = ();
    for my $analysis (@analyses) {
        for my $biotype (@biotypes) {
            push @genes,
                @{ $slice->get_all_Genes_by_type($biotype, $analysis) };
        }
    }
    return [@genes];
}

sub colour {
    my ($self, $gene, $transcript, $colours, %highlights) = @_;

    $colours ||= {};

    my $col;
    my $label;
    my $highlight;

    # Colours are set in the following order of precidence;
    # A specific colour is configured,
    # A colour is configured for the analysis logic_name
    # A colour is configured for the biotype
    # A colour is configured for the status
    my $analysis = uc($gene->analysis ? $gene->analysis->logic_name : '');
    my $biotype = uc($transcript ? $transcript->biotype : $gene->biotype);
    my $status  = uc($transcript ? $transcript->status  : $gene->status);

    my $col = ['black'];
    my $key = '';
    if (my $c = $colours->{$biotype}) {    # by biotype
        $col = $c;
        $key = $biotype;
    } elsif (my $c = $colours->{'col'}) {    # specific colour
        $col = $c;
        $key = '';
    } elsif (my $c = $colours->{$analysis}) {    # by analysis
        $col = $c;
        $key = $analysis;
    } elsif (my $c = $colours->{$status}) {      # by status
        $col = $c;
        $key = $status;
    }

    # Gramene-specific hack to highlight genes with orthologs
    if (my $homcol = $self->my_config('highlight_homologues')) {
        if (scalar(@{ $gene->get_all_homologous_Genes })) {
            if (my $c = $colours->{"${key}_HOMOLOG"}) {
                $col = $c;
            }
        }
    }

    my $stable_id = $transcript ? $transcript->stable_id : $gene->stable_id;
    my $ext_name
        = $transcript ? $transcript->external_name : $gene->external_name;

    if (exists $highlights{$stable_id}) {
        $highlight = $colours->{'superhi'};
    } elsif (exists $highlights{$ext_name}) {
        $highlight = $colours->{'superhi'};
    } elsif (exists $highlights{ $gene->stable_id() }) {
        $highlight = $colours->{'hi'};
    }

    if (ref($col) eq 'ARRAY') {
        $label = $col->[1];
        $col   = $col->[0];
    }
    $col   ||= 'black';
    $label ||= $gene->biotype;
    return ($col, $label, $highlight);
}

sub gene_colour {
    my ($self, $gene, $colours, %highlights) = @_;
    return $self->colour($gene, undef, $colours, %highlights);
}

sub href {
    my ($self, $gene, $transcript, %highlights) = @_;

    my $gid = $gene->stable_id();
    my $tid = $transcript ? $transcript->stable_id() : '';

    my $script_name = (
        $ENV{'ENSEMBL_SCRIPT'} eq 'genesnpview'
        ? 'genesnpview'
        : 'geneview'
    );

    # Check whether href is internal on gene_id or transcript_id
    if ($self->my_config('_href_only') eq '#gid'
        and exists $highlights{$gid})
    {
        return ("#$gid");
    }
    if ($self->my_config('_href_only') eq '#tid'
        and exists $highlights{$gid})
    {
        return ("#$tid");
    }

    my $species = $self->{container}{_config_file_name_};
    my $db = $self->my_config('db_alias') || 'core';
    return "/$species/$script_name?db=$db&gene=$gid";
}

sub gene_href {
    my ($self, $gene, %highlights) = @_;
    return $self->href($gene, undef, %highlights);
}

sub zmenu {
    my ($self, $gene, $transcript) = @_;

    my $sp   = $self->{container}{_config_file_name_};
    my $db   = $self->my_config('db_alias') || 'core';
    my $name = $self->my_config('db_alias') || 'Ensembl';

    my $gid   = $gene->stable_id();
    my $zmenu = {
        'caption'        => "$name Gene",
        "01:Gene:$gid"   => "/$sp/geneview?gene=$gid&db=$db",
        '04:Export Gene' => "/$sp/exportview?tab=fasta&"
            . "type=feature&ftype=gene&id=$gid",
    };

    my $top_hit = undef;
    $top_hit = $self->_get_top_protein_hit($gene);
    if (defined $top_hit) {
        $zmenu->{
            "08:Top NRAA Alignment: @{[ $top_hit->hseqname ]} (e=@{[ $top_hit->p_value ]})"
            }
            = "http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=Protein&cmd=search&term=@{[ $top_hit->hseqname ]}";
    } else {
        $zmenu->{"08:No NRAA Alignments"} = '';
    }

    if ($transcript) {
        my $tid    = $transcript->stable_id;
        my $tname  = $transcript->external_name || $tid;
        my $ext_db = $transcript->external_db || '';
        $tname = $ext_db ? "$ext_db:$tname" : $tname;
        $zmenu->{"00:$tname"}       = '';
        $zmenu->{"02:Transcr:$tid"} = "/$sp/transview?transcript=$tid&db=$db";
        $zmenu->{'05:Export cDNA'}  = "/$sp/exportview?tab=fasta&"
            . "type=feature&ftype=cdna&id=$tid";
        my $translation = $transcript->translation;

        if ($translation) {
            my $pid = $translation->stable_id;
            $zmenu->{"03:Peptide:$pid"}
                = "/$sp/protview?transcript=$tid&db=$db";
            $zmenu->{'06:Export Peptide'} = "/$sp/exportview?tab=fasta&"
                . "type=feature&ftype=peptide&id=$pid";
        }
    } else {    # No transcript
        my $gname  = $gene->external_name || $gid;
        my $ext_db = $gene->external_db   || '';
        $gname = $ext_db ? "$ext_db:$gname" : $gname;
        $zmenu->{"01:$gname"} = '';
    }

    if ($ENV{'ENSEMBL_SCRIPT'} =~ /snpview/) {
        $zmenu->{'07:Gene SNP view'} = "/$sp/genesnpview?gene=$gid&db=$db";
    }

    my $last_index = 10;
    if (my $col = $self->my_config('highlight_homologues')) {
        foreach my $homolog (@{ $gene->get_all_homologous_Genes }) {
            my ($ogene, $homology, $osp) = @$homolog;
            my $oid = $ogene->display_id;
            my $lsp = $osp;
            $lsp =~ s/\s/_/g;
            $zmenu->{ ++$last_index . ":$osp ortholog: $oid" }
                = $ogene->dbID ? "/$lsp/geneview?gene=$oid" : '';
        }
    }

    if ($gene->biotype() eq 'corrupted_translation') {
        $zmenu->{"@{[++$last_index]}:Translation inconsistent with FGENESH"}
            = '';
    }

    return $zmenu;
}

=head2 _get_top_protein_hit

    Returns the top protein hit to to the given gene

=cut

sub _get_top_protein_hit {
    my $self = shift;
    my $gene = shift;

    # Schwartzian transform
    # Get a list of peptide alignments, sorted by the p-value, minimum first
    my @peptide_alignments
        = map { $_->[1] }
        sort  { $a->[0] <=> $b->[0] }
        map { [ $_->p_value, $_ ] }
        @{ $gene->feature_Slice->get_all_ProteinAlignFeatures };
    return $peptide_alignments[0];
}

sub gene_zmenu {
    my ($self, $gene) = @_;
    return $self->zmenu($gene);
}

sub text_label {
    my ($self, $gene, $transcript) = @_;

    my $obj = $transcript || $gene || return '';

    return '' if $self->my_config('compact');

    my $tid = $obj->stable_id();
    my $eid = $obj->external_name();
    my $gid = $gene->stable_id();
    my $id  = $gid || $eid || $tid;

    my $Config = $self->{config};
    my $short_labels = $Config->get('_settings', 'opt_shortlabels');

    if ($Config->{'_both_names_'} eq 'yes') {
        $id .= $eid ? " ($eid)" : '';
    }
    if (!$Config->get('_settings', 'opt_shortlabels')) {
        my $type = (
              $gene->analysis
            ? $gene->analysis->logic_name
            : 'Generic trans.'
        );
        $id .= "\n$type";
    }
    return $id;
}

sub gene_text_label {
    my ($self, $gene) = @_;
    return $self->text_label($gene);
}

sub error_track_name { return $_[0]->my_label }

1;
