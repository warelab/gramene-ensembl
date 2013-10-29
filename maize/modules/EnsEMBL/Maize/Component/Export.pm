package EnsEMBL::Maize::Component::Export;

use EnsEMBL::Web::Component;
use EnsEMBL::Web::Form;
use EnsEMBL::Web::Document::SpreadSheet;

use Readonly;

use strict;
use warnings;

# no warnings "uninitialized";

our @ISA = qw(EnsEMBL::Web::Component);

use EnsEMBL::Maize::SeqDumper;
use CGI qw(escapeHTML);
use POSIX;

Readonly my $CONTIG_COORD_SYSTEM => 'contig';
Readonly my $CLONE_COORD_SYSTEM  => 'clone';

Readonly my @COMMON_BAC_FEATURES => (
    [ 'similarity' => 'Similarity features' ],
    [ 'repeat'     => 'Repeat features' ],
    [ 'mdr'        => 'Mathematically-Defined Repeats' ],
    [ 'improved'   => 'Improved Sequences' ],
);

Readonly my @COMMON_FPC_FEATURES => (
    [ 'clones'         => 'All FPC Clones' ],
    [ 'acc_clones'     => 'Sequenced FPC Clones' ],
    [ 'contigs'        => 'FPC Contigs' ],
    [ 'corebinmarkers' => 'Core-Bin Markers' ],
    [ 'ssr_markers'    => 'Electronic SSR Markers' ],
    [ 'overgos'        => 'Overgo Markers' ],
    [ 'other_markers'  => 'Other Markers' ],
);

Readonly my $NUCLEOTIDE_TYPE_KEY => 'nucleotide';
Readonly my $QUALITY_TYPE_KEY    => 'quality';

Readonly my $FASTA_QUALITIES_PER_LINE => 20;
Readonly my $FASTA_BASES_PER_LINE => 60;

Readonly my %FASTA_FETCHERS => (
    $NUCLEOTIDE_TYPE_KEY => \&fasta_genomic_nucleotide,
    $QUALITY_TYPE_KEY    => \&fasta_genomic_quality,
);

Readonly my @MDR_ANALYSES => map {"mdr_$_"} (0.25, 1, 2, 3);
Readonly my %MIN_COPIES => (map { $_ => ceil(10**$_) } (0.25, 1, 2, 3));

Readonly my %FORMATS => (
    'gff' => {
        'delim'  => "\t",
        'fields' => [
            qw(seqname source feature start end score strand frame attributes)
        ],
    },
    'csv' => {
        'delim'  => ',',
        'fields' => [qw(seqname refseq feature subtype start end)],
    },
    'tab' => {
        'delim'  => "\t",
        'fields' => [qw(seqname refseq feature subtype start end)],
    },
);

sub stage1 {
    my ($panel, $object) = @_;
    $panel->print(@{ [ $panel->form('stage1_form')->render() ] });
    return 1;
}

sub stage2 {
    my ($panel, $object) = @_;
    $panel->print(@{ [ $panel->form('stage2_form')->render() ] });
    return 1;
}

sub get_map_announcement_snippets {
    my ($object)      = @_;
    my $note1         = 'the FPC map';
    my $note2         = 'a specific BAC sequence map';
    my $other_species = 'Zea_mays2';
    if ($object->species ne 'Zea_mays') {
        ($note1, $note2) = ($note2, $note1);
        $other_species = 'Zea_mays';
    }
    return ($note1, $note2, $other_species);
}

sub stage1_form {
    my ($panel, $object) = @_;
    my @errors = (
        '',
        qq(Sorry, there was a problem locating the requested DNA sequence. Please check your choices - including chromosome name - and try again.),
        qq(Sorry, your chosen anchor points appear to be on different sequence regions. Please check your choices and try again.)
    );

    my ($current_map, $other_map, $other_species)
        = get_map_announcement_snippets($object);
    my $form = EnsEMBL::Web::Form->new('stage1_form',
        "/@{[$object->species]}/exportview", 'get');
    my $sitetype = ucfirst(lc($object->species_defs->ENSEMBL_SITETYPE));

    my $information = <<HTML;
<p>
Choose one or two features from the same chromosome as anchor points
and display the region between them.
Both features must be mapped to the current $sitetype golden
tiling path.  If you select "None" for the second feature,
the display will be based around the first feature.
</p>
<div class="announcement"/>
<strong>Note:</strong> You are about to export information about $current_map.
To export information about $other_map, click <a href="/$other_species/exportview">here</a>
</div>
HTML

    $form->add_element(
        'type'  => 'Information',
        'value' => $information,
    );

    if ($object->param('error')) {
        my $error_text = $errors[ $object->param('error') ];
        $form->add_element(
            'type'  => 'Information',
            'value' => '<p class="error">'
                . $error_text
                . ' If you continue to have a problem, please contact <a href="mailto:helpdesk@ensembl.org">helpdesk@ensembl.org</a>.</strong></p>'
        );
    }

    $form->add_element('type' => 'SubHeader', 'value' => 'Region');

    if (defined($object->param('l'))
        && $object->param('l') =~ /([\w\.]+):(-?[\.\w]+)-([\.\w]+)/)
    {
        my ($sr, $s, $e) = ($1, $2, $3);
        $object->input_param('seq_region_name', $sr);
        $object->input_param('type1',           'bp');
        $object->input_param('anchor1',         $s);
        $object->input_param('type2',           'bp');
        $object->input_param('anchor2',         $e);
    }
    $form->add_element(
        'type'  => 'String',
        'label' => 'Chromosome name/fragment',
        'style' => 'small',
        'name'  => 'seq_region_name',
        'id'    => 'seq_region_name',
        'value' => $object->param('seq_region_name')
    );

    my @types = (
        { 'value' => 'band',       'name' => 'Band' },
        { 'value' => 'region',     'name' => 'Region' },
        { 'value' => 'marker',     'name' => 'Marker' },
        { 'value' => 'bp',         'name' => 'Base pair' },
        { 'value' => 'gene',       'name' => 'Gene' },
        { 'value' => 'peptide',    'name' => 'Peptide' },
        { 'value' => 'transcript', 'name' => 'Transcript' },
    );
    my @types_1 = @types;
    $form->add_element(
        'type'         => 'DropDownAndString',
        'select'       => 'select',
        'name'         => 'type1',
        'label'        => 'From (type):',
        'values'       => \@types_1,
        'value'        => $object->param('type1') || 'bp',
        'string_name'  => 'anchor1',
        'string_label' => 'From (value)',
        'string_value' => $object->param('anchor1'),
        'style'        => 'medium',
        'required'     => 'yes'
    );
    unshift(@types, { 'value' => 'none', 'name' => 'None' });
    $form->add_element(
        'type'         => 'DropDownAndString',
        'select'       => 'select',
        'name'         => 'type2',
        'label'        => 'To (type):',
        'values'       => \@types,
        'value'        => $object->param('type2') || 'bp',
        'string_name'  => 'anchor2',
        'string_label' => 'To (value)',
        'style'        => 'medium',
        'string_value' => $object->param('anchor2')
    );
    $form->add_element('type' => 'SubHeader', 'value' => 'Context');
    $form->add_element(
        'type'     => 'String',
        'style'    => 'short',
        'required' => 'no',
        'value'    => '',
        'name'     => 'upstream',
        'label'    => 'Bp upstream (to the left)'
    );

    $form->add_element(
        'type'     => 'String',
        'required' => 'no',
        'value'    => '',
        'style'    => 'short',
        'name'     => 'downstream',
        'label'    => 'Bp downstream (to the right)'
    );
    $form->add_element('type' => 'SubHeader', 'value' => 'Output format');

    my $formats = $object->__data->{'formats'};
    my @formats = ();
    foreach my $super (
        sort { $formats->{$a}{'name'} cmp $formats->{$b}{'name'} }
        keys %$formats
        )
    {
        foreach my $format (
            sort {
                $formats->{$super}{'sub'}{$a}
                    cmp $formats->{$super}{'sub'}{$b}
            } keys %{ $formats->{$super}{'sub'} }
            )
        {
            push @formats,
                {
                'group' => $formats->{$super}{'name'},
                'value' => $format,
                'name'  => $formats->{$super}{'sub'}{$format}
                };
        }
    }
    $form->add_element(
        'type'   => 'DropDown',
        'select' => 'select',
        'name'   => 'format',
        'label'  => 'Output Format',
        'values' => \@formats,
        'value'  => $object->param('format')
    );
    $form->add_element(
        'type'  => 'Hidden',
        'name'  => 'action',
        'value' => 'select'
    );
    $form->add_element(
        'type'     => 'Submit',
        'value'    => 'Continue >>',
        'spanning' => 'center'
    );
    return $form;
}

sub fasta_form {
    my ($panel, $object) = @_;
    my $form = EnsEMBL::Web::Form->new('stage2_form',
        "/@{[$object->species]}/exportview", 'get');
    add_hidden_fields($form, $object);
    $form->add_element(
        'type'  => 'Hidden',
        'name'  => '_format',
        'value' => 'HTML'
    );
    my $is_transcript
        = (    $object->param('type1') eq 'transcript'
            || $object->param('type1') eq 'peptide'
            || $object->param('type1') eq 'gene')
        && ($object->param('type2') eq 'none'
        || $object->param('anchor2') eq '');

    if (!$is_transcript) {
        $form->add_element(
            'type'  => 'SubHeader',
            'value' => 'Sequence type',
        );
        $form->add_element(
            'type'  => 'DropDown',
            'name'  => 'sequence_type',
            'label' => 'FASTA type',
            'value' => $object->param('sequence_type')
                || $NUCLEOTIDE_TYPE_KEY,
            'values' => [
                { 'value' => $NUCLEOTIDE_TYPE_KEY, 'name' => 'Nucleotides' },
                { 'value' => $QUALITY_TYPE_KEY, 'name' => 'Base qualities' },
            ]
        );
    }
    $form->add_element(
        'type'  => 'SubHeader',
        'value' => 'FASTA format options'
    );

    if ($is_transcript) {
        ## We have a transcript... lets give transcript options...
        my @options = (
            [ 'genomic' => 'Genomic' ],
            [ 'cdna'    => 'cDNA' ],
            [ 'coding'  => 'Coding sequence' ],
            [ 'peptide' => 'Peptide sequence' ],
            [ 'utr5'    => "5' UTR" ],
            [ 'utr3'    => "3' UTR" ]
        );
        my %checked = map { $_ => 'yes' } $object->param('options');
        %checked = ('genomic' => 'yes') unless $object->param('options');
        $form->add_element(
            'type'   => 'MultiSelect',
            'name'   => 'options',
            'label'  => 'Sequence type to export',
            'values' => [
                map {
                    {   'value'   => $_->[0],
                        'name'    => $_->[1],
                        'checked' => $checked{ $_->[0] }
                    }
                    } @options
            ],
            'introduction' =>
                '<p>As you have requested export of a single Transcript or Translation, you can select if
                  you want to select the whole Genomic region, or just the sequence of the Transcript,
                  Translation or part thereof.</p>'
        );
    }
    add_HTML_select($form, $object);
    $form->add_element(
        'type'     => 'Submit',
        'value'    => 'Continue >>',
        'spanning' => 'center'
    );
    return $form;
}

sub flat_form {
    my ($panel, $object) = @_;
    my $form = EnsEMBL::Web::Form->new('stage2_form',
        "/@{[$object->species]}/exportview", 'get');
    add_hidden_fields($form, $object);
    $form->add_element(
        'type'  => 'Hidden',
        'name'  => '_format',
        'value' => 'HTML'
    );
    $form->add_element(
        'type'  => 'SubHeader',
        'value' => 'Flat file format options'
    );
    my @options
        = (@COMMON_BAC_FEATURES, [ 'contig' => 'Contig Information' ],);
    push @options, [ 'variation' => 'Variation features' ]
        if $object->species_defs->databases->{'ENSEMBL_VARIATION'};
    push @options, [ 'marker' => 'Marker features' ]
        if $object->species_defs->get_table_size(
        { -db => 'ENSEMBL_DB', -table => 'marker' });
    push @options, [ 'gene'     => 'Gene Information' ];
    push @options, [ 'vegagene' => 'Vega Gene Information' ]
        if $object->species_defs->databases->{'ENSEMBL_VEGA'};
    push @options, [ 'estgene' => 'EST Gene Information' ]
        if $object->species_defs->get_table_size(
        { -db => 'ENSEMBL_OTHERFEATURES', -table => 'gene' });
    my %checked
        = $object->param('options')
        ? map { $_ => 'yes' } $object->param('options')
        : ();
    $form->add_element(
        'type'   => 'MultiSelect',
        'class'  => 'radiocheck1col',
        'name'   => 'options',
        'label'  => 'Features to export',
        'values' => [
            map {
                {   'value'   => $_->[0],
                    'name'    => $_->[1],
                    'checked' => $checked{ $_->[0] }
                }
                } @options
        ]
    );
    add_HTML_select($form, $object);
    $form->add_element(
        'type'     => 'Submit',
        'value'    => 'Continue >>',
        'spanning' => 'center'
    );
    return $form;
}

sub add_hidden_fields {
    my ($form, $object) = @_;
    my @fields
        = qw(seq_region_name type1 anchor1 type2 anchor2 downstream upstream format l);
    my $flag = 1;
    foreach (@fields) {
        next unless defined $object->param($_);
        $form->add_element(
            'type'  => 'Hidden',
            'name'  => $_,
            'value' => $object->param($_)
        );
        next if $_ eq 'format';
        $flag = 0;
    }
    if ($flag) {
        $form->add_element(
            'type'  => 'Hidden',
            'name'  => 'l',
            'value' => $object->seq_region_name . ':'
                . $object->seq_region_start . '-'
                . $object->seq_region_end
        );
    }
    my $text = "<p>You are exporting @{[$object->seq_region_type_and_name]} 
                @{[ $object->thousandify( $object->seq_region_start ) ]} - 
                @{[ $object->thousandify( $object->seq_region_end ) ] }.</p>\n";
    my @O = ();
    if (   $object->param('seq_region_name')
        && $object->seq_region_name eq $object->param('seq_region_name'))
    {
        push @O, $object->seq_region_type_and_name;
    }
    push @O,
        ucfirst($object->param('type1')) . ' ' . $object->param('anchor1')
        if $object->param('anchor1');
    push @O,
        ucfirst($object->param('type2')) . ' ' . $object->param('anchor2')
        if $object->param('anchor2') && $object->param('type2') ne 'none';
    push @O, " plus " . $object->param('downstream') . " basepairs downstream"
        if $object->param('downstream');
    push @O, " plus " . $object->param('upstream') . " basepairs upstream"
        if $object->param('upstream');

    if (@O) {
        $text
            .= "<blockquote>This region is defined by: "
            . join(", ", @O)
            . "</blockquote>\n";
    }
    $form->add_element('type' => 'Information', 'value' => $text);
    $form->add_element(
        'type'  => 'Hidden',
        'name'  => 'action',
        'value' => 'export'
    );
}

sub add_HTML_select {
    my ($form, $object) = @_;
    my @options = (
        [ 'html' => 'HTML' ],
        [ 'txt'  => 'Text' ],
        [ 'gz'   => 'Compressed text (.gz)' ]
    );
    $form->add_element(
        'type'  => 'DropDown',
        'name'  => 'output',
        'label' => 'Output format',
        'value' => $object->param('output') || 'html',
        'values' =>
            [ map { { 'value' => $_->[0], 'name' => $_->[1] } } @options ]
    );
    $form->add_attribute('onSubmit',
        qq(this.elements['_format'].value='HTML';this.target='_self';flag='';for(var i=0;i<this.elements['output'].length;i++){if(this.elements['output'][i].checked){flag=this.elements['output'][i].value;}}if(flag=='txt'){this.elements['_format'].value='Text';this.target='_blank'}if(flag=='gz'){this.elements['_format'].value='TextGz';})
    );
}

sub fpc_features_form {
    my ($panel, $object) = @_;
    my $form = EnsEMBL::Web::Form->new('stage2_form',
        "/@{[$object->species]}/exportview", 'get');
    $form->add_element(
        'type'  => 'Hidden',
        'name'  => '_format',
        'value' => 'HTML'
    );
    add_hidden_fields($form, $object);
    $form->add_element(
        'type'  => 'SubHeader',
        'value' => 'Feature format options'
    );
    my @options = @COMMON_FPC_FEATURES;
    push @options, [ 'variation' => 'Variation features' ]
        if $object->species_defs->databases->{'ENSEMBL_VARIATION'};
    my %checked = map { $_ => 'yes' } $object->param('options');
    $form->add_element(
        'type'   => 'MultiSelect',
        'class'  => 'radiocheck1col',
        'name'   => 'options',
        'label'  => 'Features to export',
        'values' => [
            map {
                {   'value'   => $_->[0],
                    'name'    => $_->[1],
                    'checked' => $checked{ $_->[0] }
                }
                } @options
        ]
    );
    add_HTML_select($form, $object);
    $form->add_element(
        'type'     => 'Submit',
        'value'    => 'Continue >>',
        'spanning' => 'center'
    );
    return $form;
}

sub features_form {
    my ($panel, $object) = @_;
    my $form = EnsEMBL::Web::Form->new('stage2_form',
        "/@{[$object->species]}/exportview", 'get');
    $form->add_element(
        'type'  => 'Hidden',
        'name'  => '_format',
        'value' => 'HTML'
    );
    add_hidden_fields($form, $object);
    $form->add_element(
        'type'  => 'SubHeader',
        'value' => 'Feature format options'
    );
    my @options = @COMMON_BAC_FEATURES;
    push @options, [ 'variation' => 'Variation features' ]
        if $object->species_defs->databases->{'ENSEMBL_VARIATION'};
    push @options, [ 'gene' => 'Gene Information' ];
    my %checked = map { $_ => 'yes' } $object->param('options');
    $form->add_element(
        'type'   => 'MultiSelect',
        'class'  => 'radiocheck1col',
        'name'   => 'options',
        'label'  => 'Features to export',
        'values' => [
            map {
                {   'value'   => $_->[0],
                    'name'    => $_->[1],
                    'checked' => $checked{ $_->[0] }
                }
                } @options
        ]
    );
    add_HTML_select($form, $object);
    $form->add_element(
        'type'     => 'Submit',
        'value'    => 'Continue >>',
        'spanning' => 'center'
    );
    return $form;
}

sub features {
    my ($panel, $object) = @_;
    my @common_fields
        = qw( seqname source feature start end score strand frame );
    my %checked = map { $_ => 'yes' } $object->param('options');
    my @other_fields;
    push @other_fields, qw(hid hstart hend)
        if $checked{'similarity'} || $checked{'repeat'};
    push @other_fields, 'genscan' if $checked{'genscan'};
    push @other_fields, qw(gene_id transcript_id exon_id gene_type)
        if $checked{'gene'};
    push @other_fields, q(min_copies)
        if $checked{'mdr'};
    my %delim = ('gff' => "\t", 'csv' => ",", 'tab' => "\t");
    my $_opts = {
        'common' => \@common_fields,
        'other'  => \@other_fields,
        'delim'  => "\t",
        'format' => $object->param('format'),
        'delim'  => $delim{ $object->param('format') }
    };
    my @features = ();

    my $hid_function    = sub { $_[0]->hseqname };
    my $hstart_function = sub { $_[0]->hstart };
    my $hend_function   = sub { $_[0]->hend };
    if ($object->param('format') ne 'gff') {
        print join $delim{ $object->param('format') }, @common_fields,
            @other_fields;
        print "\n";
    }
    if ($checked{'similarity'}) {
        _add_bac_features(
            {   'target_coord_system' => $CLONE_COORD_SYSTEM,
                'type'                => 'similarity',
                'options'             => $_opts,
                'features'            => \@features,
                'get_function'        => 'get_all_SimilarityFeatures',
                'hid_function'        => $hid_function,
                'hstart_function'     => $hstart_function,
                'hend_function'       => $hend_function,
                'object'              => $object->slice,
            }
        );
    }
    if ($checked{'repeat'}) {
        _add_bac_features(
            {   'target_coord_system' => $CLONE_COORD_SYSTEM,
                'type'                => 'repeat',
                'options'             => $_opts,
                'features'            => \@features,
                'get_function'        => 'get_all_RepeatFeatures',
                'hid_function'    => sub { $_[0]->repeat_consensus->name },
                'hstart_function' => $hstart_function,
                'hend_function'   => $hend_function,
                'object'          => $object->slice,
            }
        );
    }

    # Shiran: Adding support for exporting MDRs
    if ($checked{'mdr'}) {
        my $rlevel_func = sub {
            $MIN_COPIES{ $_[0]->analysis->logic_name };
        };
        _add_bac_features(
            {   'target_coord_system' => $CLONE_COORD_SYSTEM,
                'type'                => 'mdr',
                'object'              => $object->slice,
                'options'             => $_opts,
                'features'            => \@features,
                'get_function'        => 'get_all_SimpleFeatures',
                'hid_function'        => sub { $_[0]->analysis->logic_name },
                'hstart_function'     => sub { $_[0]->start },
                'hend_function'       => sub { $_[0]->end },
                'analyses'            => [@MDR_ANALYSES],
                'extra'               => { 'min_copies' => $rlevel_func, },
            }
        );
    }

    if ($checked{'genscan'}) {
        foreach my $t (@{ $object->slice->get_all_PredictionTranscripts }) {
            foreach my $f (@{ $t->get_all_Exons }) {
                push @features,
                    _feature('pred.trans.', $_opts, $f,
                    { 'genscan' => $t->stable_id });
            }
        }
    }

    # Shiran: Adding support for improved regions
    if ($checked{'improved'}) {
        _add_bac_features(
            {   'target_coord_system' => $CLONE_COORD_SYSTEM,
                'object'              => $object->slice,
                'type'                => 'Improved Sequence',
                'options'             => $_opts,
                'features'            => \@features,
                'get_function'        => 'get_all_MiscFeatures',
                'hid_function'        => sub {q(Improved Sequence)},
                'hstart_function'     => sub { $_[0]->start },
                'hend_function'       => sub { $_[0]->end },
                'analyses'            => ['improved_sequence'],
            }
        );
    }
    if ($checked{'variation'}) {
        foreach my $f (@{ $object->slice->get_all_VariationFeatures }) {
            push @features, _feature('variation', $_opts, $f, {});
        }
    }

    if ($checked{'gene'}) {
        for my $DB (__gene_databases($object->species_defs)) {
            for my $g (@{ $object->slice->get_all_Genes(undef, $DB) }) {
                for my $t (@{ $g->get_all_Transcripts }) {
                    for my $f (map { $_->transform($CLONE_COORD_SYSTEM) }
                        @{ $t->get_all_Exons })
                    {
                        push @features,
                            +{
                            'start'   => $f->start,
                            'feature' => _feature(
                                'gene', $_opts, $f,
                                {   'exon_id'       => $f->stable_id,
                                    'transcript_id' => $t->stable_id,
                                    'gene_id'       => $g->stable_id,
                                    'gene_type'     => $g->status . '_'
                                        . $g->biotype
                                },
                                $DB eq 'vega' ? 'Vega' : 'Ensembl'
                            )
                            };
                    }
                }
            }
        }
    }
    $panel->print("<pre>");

    if ($object->param('format') ne 'gff') {
        print join $delim{ $object->param('format') }, @common_fields,
            @other_fields;
        print "\n";
    }
    for my $feature (
        map  { $_->{feature} }
        sort { $a->{start} <=> $b->{start} } @features
        )
    {
        $panel->print($feature);
    }

    $panel->print("</pre>");
}

sub bac_features {
    return features(@_);
}

sub bac_features_form {
    return features_form(@_);
}

sub fpc_features {
    my ($panel, $object) = @_;
    my %checked = map { $_ => 'yes' } $object->param('options');
    $panel->print("<pre>");

    my $format = $object->param('format');
    $panel->print(
        join($FORMATS{$format}{'delim'}, @{ $FORMATS{$format}{'fields'} }),
        "\n");

    my @features = ();
    if ($checked{'ssr_markers'}) {
        _add_features(
            {   'object'              => $object,
                'list'                => \@features,
                'function'            => 'get_all_MarkerFeatures',
                'analyses'            => [qw(ssr_marker)],
                'type'                => 'Marker',
                'target_coord_system' => 'chromosome',
            }
        );
    }
    if ($checked{'overgos'}) {
        _add_features(
            {   'object'   => $object,
                'list'     => \@features,
                'function' => 'get_all_MarkerFeatures',
                'analyses' => [qw(maize-overgos overgo-ap overgo-dupont)],
                'type'     => 'Marker',
                'target_coord_system' => 'chromosome',
            }
        );
    }
    if ($checked{'other_markers'}) {
        _add_features(
            {   'object'   => $object,
                'list'     => \@features,
                'function' => 'get_all_MarkerFeatures',
                'analyses' => [qw(Knob Telo Ribo Chloro Maize_marker)],
                'type'     => 'Marker',
                'target_coord_system' => 'chromosome',
            }
        );
    }
    if ($checked{'corebinmarkers'}) {
        _add_features(
            {   'object'              => $object,
                'list'                => \@features,
                'function'            => 'get_all_MarkerFeatures',
                'analyses'            => [qw(core_bin_marker)],
                'type'                => 'Marker',
                'target_coord_system' => 'chromosome',
            }
        );
    }
    if ($checked{'clones'}) {
        _add_features(
            {   'object'              => $object,
                'list'                => \@features,
                'function'            => 'get_all_MiscFeatures',
                'analyses'            => ['bac_map'],
                'type'                => 'Clone',
                'target_coord_system' => 'chromosome',
            }
        );
    }
    if ($checked{'acc_clones'}) {
        _add_features(
            {   'object'              => $object,
                'list'                => \@features,
                'function'            => 'get_all_MiscFeatures',
                'analyses'            => ['acc_bac_map'],
                'type'                => 'Clone',
                'target_coord_system' => 'chromosome',
            }
        );
    }
    if ($checked{'contigs'}) {
        _add_features(
            {   'object'              => $object,
                'list'                => \@features,
                'function'            => 'get_all_MiscFeatures',
                'analyses'            => ['superctgs'],
                'type'                => 'Contig',
                'target_coord_system' => 'chromosome',
            }
        );
    }
    for my $feature (sort { $a->start <=> $b->start } @features) {
        $panel->print(
            _fpc_feature(
                {   'feature' => $feature,
                    'format'  => $format,
                }
            )
        );
    }

    $panel->print("</pre>");
}

sub flat {
    my ($panel, $object) = @_;
    my $ftype = $object->seq_region_type;
    my $id    = $object->seq_region_name;
    if ($object->param('type2') eq 'none' || !$object->param('anchor2')) {
        $ftype = $object->param('type1');
        $id    = $object->param('anchor1');
    }
    my $seq_dumper = EnsEMBL::Maize::SeqDumper->new();
    my %checked = map { $_ => 'yes' } $object->param('options');
    foreach (
        qw(genscan similarity gene repeat variation contig marker mdr improved)
        )
    {
        $seq_dumper->disable_feature_type($_) unless $checked{$_};
    }
    if ($checked{'vegagene'}) {
        $seq_dumper->enable_feature_type("vegagene");
        $seq_dumper->attach_database('vega', $object->database('vega'));
    }
    if ($checked{'estgene'}) {
        $seq_dumper->enable_feature_type("estgene");
        $seq_dumper->attach_database('otherfeatures',
            $object->database('otherfeatures'));
    }
    $panel->print("<pre>");
    $seq_dumper->dump($object->slice, $object->param('format'), $panel);
    $panel->print("</pre>");
}

sub fasta_trans {
    my ($transObj, $extra) = @_;
    my $transcript = $transObj->Obj;
    my $id_type    = '';
    if ($transcript->isa('Bio::EnsEMBL::PredictionTranscript')) {
        $id_type = $transcript->analysis->logic_name;
    } elsif (0) {
        $id_type = 'externaldb';
    } else {
        $id_type = $transcript->status . '_' . $transcript->biotype;
    }
    my $out        = '';
    my %options    = map { $_ => 1 } $transObj->param('options');
    my $slice_name = '';
    my $id         = $transcript->stable_id;

    $id = "$extra:$id" if $extra;
    $out .= format_fasta("$id cdna:$id_type $slice_name",
        $transcript->spliced_seq)
        if $options{'cdna'};
    $out .= format_fasta("$id cds:$id_type $slice_name",
        $transcript->translateable_seq)
        if $options{'coding'} && $transcript->translation;
    if ($options{'peptide'} && $transcript->translation) {
        my $pep_id = $transcript->translation->stable_id;
        $out .= format_fasta("$id peptide:$pep_id pep:$id_type $slice_name",
            $transcript->translate->seq);
    }
    $out .= format_fasta("$id utr3:$id_type $slice_name",
        $transcript->three_prime_utr->seq)
        if $options{'utr3'} && $transcript->three_prime_utr;
    $out .= format_fasta("$id utr5:$id_type $slice_name",
        $transcript->five_prime_utr->seq)
        if $options{'utr5'} && $transcript->five_prime_utr;
    return $out;
}

sub format_fasta {
    my ($line1, $seq) = @_;
    $seq =~ s/(.{$FASTA_BASES_PER_LINE})/$1\n/g;
    return ">$line1\n$seq\n";
}

sub format_quality {
    my ($description, $qualities) = @_;
    my @output = (">$description");
    while ($qualities =~ m/((?:\d+\s?){1,$FASTA_QUALITIES_PER_LINE})/g) {
        push @output, join('', map { sprintf('%3d', $_) } split(' ', $1));
    }
    return join("\n", @output);
    
}

sub fasta {
    my ($panel, $object) = @_;
    my $seq;
    my $id   = $object->seq_region_name;
    my $desc = "dna:@{[$object->seq_region_type]} @{[$object->slice->name]}";
    ## First of all what sort of object do we have?
    my $genomic = 1;
    my $output  = '';
    foreach my $transObj (@{ $object->__data->{'transcript'} || [] }) {
        $output .= fasta_trans($transObj);
        $genomic = 0;
    }
    foreach my $geneObj (@{ $object->__data->{'gene'} || [] }) {
        foreach my $transObj (@{ $geneObj->get_all_transcripts }) {
            $output .= fasta_trans($transObj, $geneObj->stable_id);
            $genomic = 0;
        }
    }
    my %options
        = map { $_ => 1 } grep { defined($_) } $object->param('options');
    if ($genomic || $options{'genomic'}) {
        my $sequence_type = $object->param('sequence_type');
        my $fetcher       = $FASTA_FETCHERS{$sequence_type};
        $fetcher ||= \&fasta_genomic_nucleotide;
        $output .= $fetcher->($object);
    }
    $panel->print("<pre>$output</pre>");
}

=head2 fasta_genomic_nucleotide

    Returns the nucleotide FASTA sequence for a genomic object (as opposed
    to gene or transcript object)

=cut

sub fasta_genomic_nucleotide {
    my ($object) = @_;
    return format_fasta(
        "@{[$object->seq_region_name]} dna:@{[$object->seq_region_type]} @{[$object->slice->name]}",
        $object->slice->seq
    );
}

use Bio::SeqIO;
use Bio::Seq::PrimaryQual;
use IO::String;

=head2 fasta_genomic_quality

    Returns the quality FASTA sequence for a genomic object

=cut

sub fasta_genomic_quality {
    my ($object)    = @_;
    my $slice       = $object->slice;
    my $description = join(q{},
        $object->seq_region_name, ' Phred Quality ',
        $slice->name);

    my $quality_score_adaptor
        = $slice->adaptor->db->get_adaptor('QualityScore');
    my $quality_score
        = $quality_score_adaptor->fetch_raw_scores_for_Slice($slice);

    my $output = format_quality($description, $quality_score->scores());
    return $output;
}

sub pip_form {
    my ($panel, $object) = @_;
    my $form = EnsEMBL::Web::Form->new('stage2_form',
        "/@{[$object->species]}/exportview", 'get');
    my $FN = $object->temp_file_name(undef, 'XXX/X/X/XXXXXXXXXXXXXXX');
    my $seq_file = $object->species_defs->ENSEMBL_TMP_DIR_IMG . "/$FN.fa";
    my ($PATH, $FILE)
        = $object->make_directory(
        $object->species_defs->ENSEMBL_TMP_DIR_IMG . "/$FN");
    my $anno_file = $object->species_defs->ENSEMBL_TMP_DIR_IMG . "/$FN.txt";
    my $arc_file = $object->species_defs->ENSEMBL_TMP_DIR_IMG . "/$FN.tar.gz";
    my $seq_url  = $object->species_defs->ENSEMBL_TMP_URL_IMG . "/$FN.fa";
    my $anno_url = $object->species_defs->ENSEMBL_TMP_URL_IMG . "/$FN.txt";
    my $arc_url  = $object->species_defs->ENSEMBL_TMP_URL_IMG . "/$FN.tar.gz";
    open O, ">$seq_file";
    print O "T";
    close O;
    pip_seq_file($seq_file, $object);
    pip_anno_file($anno_file, $object);
    system("cd $PATH; tar cf - $FILE.fa $FILE.txt | gzip -9 > $FILE.tar.gz");
    $form->add_element(
        'type'  => 'SubHeader',
        'value' => 'PIP format options'
    );
    $form->add_element(
        'type'  => 'Information',
        'value' => qq(<p>
    Your export has been processed successfully. Please download
    the exported data by following the links below.
  </p>
  <ul>
    <li><strong>Sequence data:</strong> <a target="_blank" href="$seq_url">$FILE.fa</a> [FASTA format]</li>
    <li><strong>Annotation data:</strong> <a target="_blank" href="$anno_url">$FILE.txt</a> [@{[$object->param('format')]} format]</li>
  </ul>
  <p><strong>OR</strong>
  <ul>
    <li><strong>Combined file:</strong> <a href="$arc_url">$FILE.tar.gz</a></li>
  </ul>
  )
    );
    return $form;
}

sub pip_seq_file {
    my ($FN, $object) = @_;
    open O, ">$FN";
    (my $seq = $object->slice->seq) =~ s/(.{60})/$1\n/g;
    print O ">@{[$object->slice->name]}\n$seq";
    close O;
}

sub pip_anno_file {
    my ($FN, $object) = @_;
    open O, ">$FN";
    my $format = $object->param('format');    # either pipmaker or zpicture!
    my $slice  = $object->slice;
    my $slice_length = $slice->length;
    my $exonfunct    = "pip_anno_file_$format";
    foreach my $gene (@{ $slice->get_all_Genes(undef, undef, 1) || [] }) {

        # only include genes that don't overlap slice boundaries
        next if ($gene->start < 1 or $gene->end > $slice_length);
        my $gene_header = join(" ",
            ($gene->strand == 1 ? ">" : "<"),
            $gene->start, $gene->end,
            $gene->external_name || $gene->stable_id);
        $gene_header .= "\n";
        foreach my $transcript (@{ $gene->get_all_Transcripts }) {

            # get UTR/exon lines
            my @exons = @{ $transcript->get_all_Exons };
            @exons = reverse @exons if ($gene->strand == -1);
            no strict 'refs';
            my $out = &$exonfunct($transcript, \@exons);

            # write output to file if there are exons in the exported region
            if ($out) {
                print O $gene_header, $out;
            }
        }
    }
    close O;
}

sub pip_anno_file_vista {
    my ($transcript, $exons) = @_;
    my $coding_start = $transcript->coding_region_start;
    my $coding_end   = $transcript->coding_region_end;
    my $out;
    foreach my $exon (@{$exons}) {
        if (!$coding_start) {    # no coding region at all
            $out .= join(" ", $exon->start, $exon->end, "UTR\n");
        } elsif ($exon->start < $coding_start) {    # we begin with an UTR
            if ($coding_start < $exon->end)
            {    # coding region begins in this exon
                $out .= join(" ", $exon->start,  $coding_start - 1, "UTR\n");
                $out .= join(" ", $coding_start, $exon->end,        "exon\n");
            } else {    # UTR until end of exon
                $out .= join(" ", $exon->start, $exon->end, "UTR\n");
            }
        } elsif ($coding_end < $exon->end) {    # we begin with an exon
            if ($exon->start < $coding_end)
            {    # coding region ends in this exon
                $out .= join(" ", $exon->start,    $coding_end, "exon\n");
                $out .= join(" ", $coding_end + 1, $exon->end,  "UTR\n");
            } else {    # UTR (coding region has ended in previous exon)
                $out .= join(" ", $exon->start, $exon->end, "UTR\n");
            }
        } else {    # coding exon
            $out .= join(" ", $exon->start, $exon->end, "exon\n");
        }
    }
    return $out;
}

sub pip_anno_file_pipmaker {
    my ($transcript, $exons) = @_;
    my $coding_start = $transcript->coding_region_start;
    my $coding_end   = $transcript->coding_region_end;
    my $out;

    # do nothing for non-coding transcripts
    return unless ($coding_start);

    # add UTR line
    if (   $transcript->start < $coding_start
        or $transcript->end > $coding_end)
    {
        $out .= join(" ", "+", $coding_start, $coding_end, "\n");
    }

    # add exon lines
    foreach my $exon (@{$exons}) {
        $out .= join(" ", $exon->start, $exon->end, "\n");
    }
    return $out;
}

=pod

=head2 _add_features
    Add features

=cut

sub _add_features {
    my ($params)            = @_;
    my $object              = $params->{'object'};
    my $function            = $params->{'function'};
    my @analyses            = @{ $params->{'analyses'} };
    my $type                = $params->{'type'};
    my $target_coord_system = $params->{'target_coord_system'};
    my $features            = $params->{'list'};
    if (scalar @analyses == 0) {
        push @analyses, undef;
    }
    for my $analysis (@analyses) {

        # Transform function places the feature on the entire chromosome
        for my $feature (map { $_->transform($target_coord_system) }
            @{ $object->slice->$function($analysis) })
        {
            $feature->{__feature_type} = $type;
            if (defined $feature->analysis) {
                $feature->{__feature_subtype}
                    = $feature->analysis->logic_name;
            } elsif (defined $feature->get_scalar_attribute('embl_acc')) {
                $feature->{__accession}
                    = $feature->get_scalar_attribute('embl_acc');
                $feature->{__seqstatus}
                    = $feature->get_scalar_attribute('external') eq 'true'
                    ? 'EXTERNAL'
                    : $feature->get_scalar_attribute('seqstatus');
            }
            push @$features, $feature;
        }
    }
}

sub _fpc_feature {
    my ($params) = @_;
    my $feature  = $params->{'feature'};
    my $format   = $params->{'format'};

    my $feature_attributes = +{
        'seqname' => $feature->display_id || 'SEQ',
        'feature' => $feature->{__feature_type},
        'subtype' => $feature->{__feature_subtype},
        'start'   => $feature->start,
        'end'     => $feature->end,
        'source'  => 'Ensembl',
        'strand'  => '+',
        'frame'   => '.',
        'score'   => '.',
        'refseq'  => $feature->slice->desc,
        'attributes' => 'RefSeq "' . $feature->slice->desc . '"',
        'status'     => $feature->{__seqstatus},
        'accession'  => $feature->{__accession},
    };

    if ($format eq 'gff') {
        $feature_attributes->{'attributes'} = join('; ',
            map      {"$_=\"@{[$feature_attributes->{$_}]}\""}
                grep { defined $feature_attributes->{$_} }
                qw(refseq subtype accession status));
    }

    return (
        join(
            $FORMATS{$format}{'delim'},
            map { $feature_attributes->{$_} } @{ $FORMATS{$format}{'fields'} }
        ),
        "\n"
    );
}

sub _add_bac_features {
    my ($params)            = @_;
    my $type                = $params->{'type'};
    my $target_coord_system = $params->{'target_coord_system'};
    my $analyses_ref        = $params->{'analyses'};
    my $get_function        = $params->{'get_function'};
    my $object              = $params->{'object'};
    my $hid_function        = $params->{'hid_function'};
    my $hstart_function     = $params->{'hstart_function'};
    my $hend_function       = $params->{'hend_function'};
    my $features            = $params->{'features'};
    my $options             = $params->{'options'};
    my $extra               = $params->{'extra'};

    my @analyses = $analyses_ref ? @$analyses_ref : ();
    if (scalar @analyses == 0) {
        push @analyses, undef;
    }
    for my $analysis (@analyses) {

        # Transform function places the feature on the entire chromosome
        for my $feature (map { $_->transform($target_coord_system) }
            @{ $object->$get_function($analysis) })
        {
            my %extra = ();
            for my $parameter (keys %$extra) {
                $extra{$parameter} = &{ $extra->{$parameter} }($feature);
            }
            my $start = $feature->can('start') ? $feature->start : 1;
            my $export_feature = _feature(
                $type, $options, $feature,
                {   'hid'    => &$hid_function($feature),
                    'hstart' => &$hstart_function($feature),
                    'hend'   => &$hend_function($feature),
                    %extra,
                }
            );

            push @$features,
                +{ 'start' => $start, 'feature' => $export_feature };
        }
    }

}

sub _feature {
    my ($type, $options, $feature, $extra, $def_source) = @_;
    my $score = $feature->can('score') ? $feature->score : undef;
    $score ||= '.';
    my $frame = $feature->can('frame') ? $feature->frame : undef;
    $frame ||= '.';
    my $strand = $feature->can('strand') ? $feature->strand : undef;
    $strand ||= '.';
    $strand = '+' if $strand eq 1;
    $strand = '-' if $strand eq -1;
    my $name
        = $feature->can('entire_seq')
        && $feature->entire_seq
        ? $feature->entire_seq->name
        : undef;
    $name = $feature->seqname
        if !$name && $feature->can('seqname');
    $name ||= 'SEQ';
    $name =~ s/\s/_/g;
    my $source
        = $feature->can('source_tag')
        ? $feature->source_tag
        : undef;
    $source ||= $def_source || 'Ensembl';
    $source =~ s/\s/_/g;
    my $tag
        = $feature->can('primary_tag')
        ? $feature->primary_tag
        : undef;
    $tag ||= ucfirst(lc($type)) || '.';
    $tag =~ s/\s/_/g;
    my $start = $feature->can('start') ? $feature->start : undef;
    my $end   = $feature->can('end')   ? $feature->end   : undef;

    my @results
        = ($name, $source, $tag, $start, $end, $score, $strand, $frame);

    if ($options->{'format'} eq 'gff') {
        push @results, join "; ",
            map { defined $extra->{$_} ? "$_=$extra->{$_}" : () }
            @{ $options->{'other'} };
    } else {
        push @results, map { $extra->{$_} } @{ $options->{'other'} };
    }
    return join($options->{'delim'}, @results) . "\n";
}

sub __gene_databases {
    my $species_defs = shift;
    my @return       = ('core');
    push @return, 'vega'
        if $species_defs->databases->{'ENSEMBL_VEGA'};
    push @return, 'est'
        if $species_defs->databases->{'ENSEMBL_OTHERFEATURES'};
    return @return;
}

1;
