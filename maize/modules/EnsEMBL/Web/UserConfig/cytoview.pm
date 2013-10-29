package EnsEMBL::Web::UserConfig::cytoview;
use strict;
use EnsEMBL::Web::UserConfig;
use vars qw(@ISA);
@ISA = qw(EnsEMBL::Web::UserConfig);

sub init {
    my ($self) = @_;
    $self->{'_userdatatype_ID'} = 20;
    $self->{'_das_offset'}      = '5080';

    $self->{'general'}->{'cytoview'} = {
        '_artefacts' => [
            qw(
                assemblyexception
                bacends
                bacs
                bac_bands
                bac_map
                blast_new
                cloneset_1mb
                cloneset_32k
                cloneset_37k
                encode
                clone_legend
                gap
                gene_legend
                haplotype_links
                missing
                nod_bacs
                ntcontigs
                qtl
                ruler
                scalebar
                supercontigs
                tilepath
                cyto_bac_genes
                corebins
                )
        ],
        '_settings' => {
            'URL'              => '',
            'width'            => 900,
            'default_vc_size'  => 5000000,
            'band_box'         => 'show',
            'show_cytoview'    => 'yes',
            'imagemap'         => 1,
            'opt_pdf'          => 0,
            'opt_svg'          => 0,
            'opt_postscript'   => 0,
            'opt_lines'        => 1,
            'opt_empty_tracks' => 0,
            'opt_gene_labels'  => 1,
            'opt_zmenus'       => 1,
            'opt_zclick'       => 1,
            'bgcolor'          => 'background1',
            'bgcolour1'        => 'background2',
            'bgcolour2'        => 'background3',
            'show_bands_nav'   => 'yes',
            'zoom_gifs'        => {
                zoom1 => 200000,
                zoom2 => 500000,
                zoom3 => 1000000,
                zoom4 => 2000000,
                zoom5 => 5000000,
                zoom6 => 10000000,
                zoom7 => 20000000,
                zoom8 => 50000000
            },
            'navigation_options' =>
                [ '2mb', '1mb', 'window', 'half', 'zoom' ],
            'compara'  => [],
            'features' => [
                [ 'tilepath'       => 'Sequenced BACs' ],
                [ 'bac_map'        => 'BAC map' ],
                [ 'supercontigs'   => 'FPC Contigs' ],
                [ 'Marker'         => 'Hybridized Markers' ],
                [ 'qtl'            => 'QTLs' ],
                [ 'corebins'       => 'Virtual Core Bins' ],
                [ 'cyto_bac_genes' => 'BAC Genes' ],
            ],
            'options' => [
                [ 'nod_bacs'         => 'Nod BACs' ],
                [ 'bac_bands'        => 'Band BACs' ],
                [ 'haplotype_links'  => 'Haplotype blocks' ],
                [ 'bacs'             => 'BACs' ],
                [ 'bacends'          => 'BACends' ],
                [ 'ntcontigs'        => 'NT Contigs' ],
                [ 'ruler'            => 'Ruler' ],
                [ 'scalebar'         => 'Scale Bar' ],
                [ 'cloneset_1mb'     => '1Mb Cloneset' ],
                [ 'cloneset_37k'     => '37k Cloneset' ],
                [ 'cloneset_32k'     => '32k Cloneset' ],
                [ 'encode'           => 'Encode regions' ],
                [ 'ssr'              => 'SSRs' ],
                [ 'opt_lines'        => 'Show register lines' ],
                [ 'opt_empty_tracks' => 'Show empty tracks' ],
                [ 'opt_zmenus'       => 'Show popup menus' ],
                [ 'opt_zclick'       => '... popup on click' ],
                [ 'gap'              => 'Gaps' ],
            ],
            'menus' =>
                [qw(features DAS options export jumpto resize)]    # repeats
        },
        'stranded_contig' => {
            'on'                   => "on",
            'pos'                  => '0',
            'col'                  => 'black',
            'threshold_navigation' => '10000'
        },

## Blast and SSAHA tracks displayed if linked to from Blast/SSAHA...
## These get put beside the central track and so have low pos

        'redbox' => {
            'on'     => 'off',
            'pos'    => '1000000',
            'col'    => 'gold',
            'zindex' => -20,
        },

        'blast_new' => {
            'on'          => "on",
            'pos'         => '7',
            'col'         => 'red',
            'dep'         => '6',
            'str'         => 'b',
            'force_cigar' => 'yes',
        },

        'cloneset_1mb' => {
            'on'      => "on",
            'pos'     => '1005',
            'colours' => {
                'col_CES_AVC_MISMATCH'                        => 'red',
                'col_CES_ONLY'                                => 'grey50',
                'col_CES_UNVERIFIED'                          => 'grey50',
                'col_CLONE_ACCESSION'                         => 'gold',
                'col_CLONE_ACCESSION_END_SEQ_UCSC'            => 'gold',
                'col_END_SEQ_UCSC'                            => 'gold',
                'col_REPICK_CLONE_ACCESSION'                  => 'orange',
                'col_REPICK_CLONE_ACCESSION'                  => 'orange',
                'col_REPICK_END_CLONE_ACCESSION'              => 'orange',
                'col_REPICK_END_CLONE_ACCESSION_END_SEQ_UCSC' => 'orange',
                'col_REPICK_END_ONLY_CLONE_ACCESSION'         => 'orange',
                'col_REPICK_END_ONLY_CLONE_ACCESSION_END_SEQ_UCSC' =>
                    'orange',
                'col_REPICK_END_SEQ_UCSC'          => 'orange',
                'col_SSAHA2'                       => 'contigblue2',
                'col_TELOMERE'                     => 'grey50',
                'col_TPF_CLONE_ACCESSION'          => 'gold',
                'col_TPF_CLONE_ACCESSION_2'        => 'gold',
                'lab_CES_AVC_MISMATCH'             => 'white',
                'lab_CES_ONLY'                     => 'black',
                'lab_CES_UNVERIFIED'               => 'black',
                'lab_CLONE_ACCESSION'              => 'black',
                'lab_CLONE_ACCESSION_END_SEQ_UCSC' => 'black',
                'lab_END_SEQ_UCSC'                 => 'black',
                'lab_REPICK_CLONE_ACCESSION'       => 'white',
                'lab_REPICK_CLONE_ACCESSION'       => 'white',
                'lab_REPICK_END_CLONE_ACCESSION'   => 'white',
                'lab_REPICK_END_CLONE_ACCESSION_END_SEQ_UCSC'      => 'white',
                'lab_REPICK_END_ONLY_CLONE_ACCESSION'              => 'white',
                'lab_REPICK_END_ONLY_CLONE_ACCESSION_END_SEQ_UCSC' => 'white',
                'lab_REPICK_END_SEQ_UCSC'                          => 'white',
                'lab_SSAHA2'                                       => 'black',
                'lab_TELOMERE'                                     => 'black',
                'lab_TPF_CLONE_ACCESSION'                          => 'black',
                'lab_TPF_CLONE_ACCESSION_2'                        => 'black',
                'seq_len'                                          => 'black',
                'fish_tag'                                         => 'black',
            },
            'str'                  => 'r',
            'dep'                  => '9999',
            'threshold_navigation' => '10000000',
            'fish'                 => 'FISH',
            'available'            => 'features mapset_cloneset_1mb',
        },
        'cloneset_37k' => {
            'on'      => "on",
            'pos'     => '1006',
            'colours' => {
                'col_CES_AVC_MISMATCH'                        => 'red',
                'col_CES_ONLY'                                => 'grey50',
                'col_CES_UNVERIFIED'                          => 'grey50',
                'col_CLONE_ACCESSION'                         => 'gold',
                'col_CLONE_ACCESSION_END_SEQ_UCSC'            => 'gold',
                'col_END_SEQ_UCSC'                            => 'gold',
                'col_REPICK_CLONE_ACCESSION'                  => 'orange',
                'col_REPICK_CLONE_ACCESSION'                  => 'orange',
                'col_REPICK_END_CLONE_ACCESSION'              => 'orange',
                'col_REPICK_END_CLONE_ACCESSION_END_SEQ_UCSC' => 'orange',
                'col_REPICK_END_ONLY_CLONE_ACCESSION'         => 'orange',
                'col_REPICK_END_ONLY_CLONE_ACCESSION_END_SEQ_UCSC' =>
                    'orange',
                'col_REPICK_END_SEQ_UCSC'          => 'orange',
                'col_SSAHA2'                       => 'contigblue2',
                'col_TELOMERE'                     => 'grey50',
                'col_TPF_CLONE_ACCESSION'          => 'gold',
                'col_TPF_CLONE_ACCESSION_2'        => 'gold',
                'lab_CES_AVC_MISMATCH'             => 'white',
                'lab_CES_ONLY'                     => 'black',
                'lab_CES_UNVERIFIED'               => 'black',
                'lab_CLONE_ACCESSION'              => 'black',
                'lab_CLONE_ACCESSION_END_SEQ_UCSC' => 'black',
                'lab_END_SEQ_UCSC'                 => 'black',
                'lab_REPICK_CLONE_ACCESSION'       => 'white',
                'lab_REPICK_CLONE_ACCESSION'       => 'white',
                'lab_REPICK_END_CLONE_ACCESSION'   => 'white',
                'lab_REPICK_END_CLONE_ACCESSION_END_SEQ_UCSC'      => 'white',
                'lab_REPICK_END_ONLY_CLONE_ACCESSION'              => 'white',
                'lab_REPICK_END_ONLY_CLONE_ACCESSION_END_SEQ_UCSC' => 'white',
                'lab_REPICK_END_SEQ_UCSC'                          => 'white',
                'lab_SSAHA2'                                       => 'black',
                'lab_TELOMERE'                                     => 'black',
                'lab_TPF_CLONE_ACCESSION'                          => 'black',
                'lab_TPF_CLONE_ACCESSION_2'                        => 'black',
                'seq_len'                                          => 'black',
                'fish_tag'                                         => 'black',
            },
            'str'                  => 'r',
            'dep'                  => '9999',
            'threshold_navigation' => '10000000',
            'fish'                 => 'FISH',
            'available'            => 'features mapset_cloneset_37k',
        },

        'cloneset_32k' => {
            'on'                   => 'on',
            'pos'                  => '1007',
            'colour'               => 'green',
            'str'                  => 'r',
            'dep'                  => '9999',
            'threshold_navigation' => '10000000',
            'available'            => 'features mapset_cloneset_32k'
        },

        'encode' => {
            'on'                   => 'on',
            'pos'                  => '1010',
            'colour'               => 'salmon',
            'label'                => 'black',
            'str'                  => 'r',
            'dep'                  => '9999',
            'threshold_navigation' => '10000000',
            'available'            => 'features mapset_encode'
        },

        'haplotype_links' => {
            'on'                   => "on",
            'pos'                  => '999',
            'col'                  => 'red',
            'lab'                  => 'white',
            'available'            => 'features mapset_haplotype',
            'str'                  => 'r',
            'dep'                  => '9999999',
            'threshold_navigation' => '10000000',
            'outline_threshold'    => '35000000'
        },
        'nod_bacs' => {
            'on'                   => "on",
            'pos'                  => '997',
            'col'                  => 'red',
            'lab'                  => 'black',
            'available'            => 'features mapset_nod_bacs',
            'str'                  => 'r',
            'dep'                  => '9999999',
            'threshold_navigation' => '100000',
            'outline_threshold'    => '350000'
        },
        'bac_bands' => {
            'on'        => "on",
            'pos'       => '996',
            'col'       => 'darkred',
            'lab'       => 'black',
            'available' => 'features mapset_bacs_bands',
            'colours'   => {
                'col_unmapped'   => 'contigblue2',
                'col_conflict'   => 'darkslateblue',
                'col_consistent' => 'springgreen4',
                'lab_unmapped'   => 'white',
                'lab_conflict'   => 'white',
                'lab_consistent' => 'white'
            },
            'dep'               => '9999',
            'str'               => 'r',
            'outline_threshold' => '350000'
        },

        'bac_map' => {
            'on'        => 'off',
            'pos'       => '995',
            'col'       => 'green',
            'font'      => 'Tiny',
            'lab'       => 'black',
            'available' => 'features mapset_bac_map',
            'colours'   => {
                'col_Free'        => 'grey70',        #'gray80',
                'col_Phase0Ac'    => 'thistle2',
                'col_Committed'   => 'mediumpurple1',
                'col_PreDraftAc'  => 'plum',
                'col_Redundant'   => 'gray80',
                'col_Reserved'    => 'gray80',
                'col_DraftAc'     => 'gold2',
                'col_FinishAc'    => 'gold3',
                'col_Abandoned'   => 'gray80',
                'col_FULLTOP'     => 'lightskyblue',
                'col_PREFIN'      => 'royalblue1',
                'col_ACTIVEFIN'   => 'royalblue3',
                'col_IMPROVED'    => 'royalblue4',
                'col_External'    => 'peru',
                'col_Accessioned' => 'red',           # Maize. Ori='thistle2',
                'col_Unknown'     => 'gray80',
                'col_'            => 'gray80',
                'lab_Free'        => 'black',
                'lab_Phase0Ac'    => 'black',
                'lab_Committed'   => 'black',
                'lab_PreDraftAc'  => 'black',
                'lab_Redundant'   => 'black',
                'lab_Reserved'    => 'black',
                'lab_DraftAc'     => 'black',
                'lab_FinishAc'    => 'black',
                'lab_Abandoned'   => 'black',
                'lab_Accessioned' => 'black',
                'lab_Unknown'     => 'black',
                'lab_'            => 'black',
                'bacend'          => 'black',
                'seq_len'         => 'black',
            },
            'str'                  => 'r',
            'dep'                  => '9999999',
            'threshold_navigation' => '100000',
            'full_threshold'       => '50000',
            'outline_threshold'    => '350000'
        },

        'supercontigs' => {
            'on'        => "on",
            'pos'       => '990',
            'col'       => 'green',
            'lab'       => 'black',
            'available' => 'features mapset_superctgs',
            'colours'   => {
                'col1' => 'darkgreen',
                'col2' => 'green',
                'lab1' => 'white',
                'lab2' => 'black',
            },
            'str'                  => 'r',
            'dep'                  => '9999999',
            'threshold_navigation' => '10000000'
        },

        'ntcontigs' => {
            'on'        => "off",
            'pos'       => '991',
            'col'       => 'green',
            'lab'       => 'black',
            'available' => 'features mapset_ntctgs',
            'colours'   => {
                'col1' => 'darkgreen',
                'col2' => 'green',
                'lab1' => 'black',
                'lab2' => 'black',
            },
            'str'                  => 'r',
            'dep'                  => '0',
            'threshold_navigation' => '10000000'
        },

        'corebins' => {
            'on'        => 'on',
            'pos'       => '500',
            'col'       => 'purple',
            'lab'       => 'black',
            'available' => 'features mapset_core_bins',
            'str'       => 'r',
            'dep'       => '9999999'
        },

        'tilepath' => {
            'on'        => "on",
            'pos'       => '1011',
            'col'       => 'green',
            'lab'       => 'black',
            'available' => 'features mapset_acc_bac_map',
            'colours'   => {
                'col_FULLTOP'   => 'lightskyblue',
                'col_PREFIN'    => 'royalblue1',
                'col_ACTIVEFIN' => 'royalblue3',
                'col_IMPROVED'  => 'royalblue4',
                'col_External'  => 'peru',
                'label'         => 'black',
                'bacend'        => 'black',
                'seq_len'       => 'black',
            },
            'fish'                 => 'FISH',
            'str'                  => 'r',
            'dep'                  => '9999999',
            'threshold_navigation' => '10000000',
            'outline_threshold'    => '350000'
        },
        'marker' => {
            'on'        => "on",
            'pos'       => '1501',
            'dep'       => '20',
            'str'       => 'r',
            'colours'   => { $self->{'_colourmap'}->colourSet('marker') },
            'labels'    => 'on',
            'available' => 'features marker',
        },
        'marker_overgo' => {
            'on'     => "on",
            'pos'    => '2000',
            'dep'    => '200',
            'str'    => 'r',
            'col'    => 'green',
            'labels' => 'on',
            'available' =>
                'features OVERGO',    ## track will work with or without
        },

        'gap' => {
            'on'        => "off",
            'pos'       => '8020',
            'col1'      => 'red',
            'col2'      => 'orange',
            'lab1'      => 'black',
            'lab2'      => 'black',
            'available' => 'features mapset_gap',
            'str'       => 'r',
        },
        'qtl' => {
            'on'        => 'on',
            'pos'       => '1504',
            'col'       => 'lightcoral',
            'lab'       => 'black',
            'available' => 'features qtl',
            'dep'       => '99999',
            'str'       => 'r',
        },

        'scalebar' => {
            'on'         => "on",
            'nav'        => "on",
            'pos'        => '8000',
            'col'        => 'black',
            'str'        => 'b',
            'abbrev'     => 'on',
            'navigation' => 'on'
        },

        'sub_repeat' => {
            'on'                   => "on",
            'pos'                  => '4087',
            'str'                  => 'r',
            'col'                  => 'gray50',
            'threshold'            => '2000',
            'navigation_threshold' => '1000',
            'navigation'           => 'on'
        },
        'ruler' => {
            'on'  => "on",
            'pos' => '9010',
            'col' => 'black',
        },
        'gene_legend' => {
            'on'  => "on",
            'str' => 'r',
            'pos' => '100000',
            'src' => 'all',      # 'ens' or 'all'
            'dep' => '6',
        },
        'clone_legend' => {
            'on'  => "on",
            'str' => 'r',
            'pos' => '100000',
            'src' => 'all',      # 'ens' or 'all'
            'dep' => '6',
        },
        'missing' => {
            'on'  => "on",
            'str' => 'r',
            'pos' => '100001',
            'src' => 'all',      # 'ens' or 'all'
            'dep' => '6',
        },
        'bacends' => {
            'on'        => "off",
            'pos'       => '1025',
            'col'       => 'red',
            'lab'       => 'black',
            'available' => 'features bacends',
            'dep'       => 6,
            'str'       => 'r'
        },
        'urlfeature' => {
            'on'                   => "on",
            'pos'                  => '7099',
            'str'                  => 'b',
            'col'                  => 'red',
            'force_cigar'          => 'yes',
            'dep'                  => 6,
            'navigation'           => 'on',
            'navigation_threshold' => '2000',
            'threshold'            => '2000',
        },

        'bacs' => {
            'on'        => "off",
            'pos'       => '1020',
            'col'       => 'green1',
            'lab'       => 'black',
            'available' => 'features mapset_bacs',
            'colours'   => {
                'col_unmapped'   => 'cadetblue3',
                'col_conflict'   => 'firebrick1',
                'col_consistent' => 'darkgreen',
                'lab_unmapped'   => 'white',
                'lab_conflict'   => 'white',
                'lab_consistent' => 'white'
            },
            'dep'               => '9999',
            'str'               => 'r',
            'outline_threshold' => '350000'
        },
        'assemblyexception' => {
            'on'         => "on",
            'pos'        => '9998',
            'str'        => 'x',
            'lab'        => 'black',
            'dep'        => 6,
            'navigation' => 'on',
        },

    };
    $self->ADD_GENE_TRACKS();
    $self->ADD_SYNTENY_TRACKS();
    $self->ADD_ALL_MARKER_FEATURES();
    $self->ADD_DNA_FEATURES();
    $self->add_legend(
        'track'   => 'tilepath',
        'name'    => 'clone_legend',
        'caption' => 'BAC Legend',
        'columns' => 4,
    );
    $self->add_legend(
        'track'   => 'Markers',
        'name'    => 'marker_legend',
        'caption' => 'Marker Legend',
        'columns' => 4,
    );
}

=head2 ADD_DNA_FEATURES

    Add DNA alignments

=cut

sub ADD_DNA_FEATURES {
    my $self                   = shift;
    my $POS                    = shift || 1013;
    my %dna_feature_parameters = (
        'genes' => {
            'subtype_labeler' => sub {
                $self->biotype_label(@_);
            },
            'grouper' => sub {
                [   qw/protein_coding protein_coding_unsupported
                        transposon_pseudogene/
                ];
            },
        },
        'alignments' => {
            'subtype_labeler' => sub {
                $self->logic_name_label(@_);
            },
            'grouper' => sub {
                $self->grouping(@_);
            },
        },
        'repeats' => {
            'subtype_labeler' => sub {
                my ($label) = ($_[0] =~ m|^MIPS/REcat (.*)|);
                return $label;
            },
            'grouper' => sub {
                [   'MIPS/REcat Class I Retroelements',
                    'MIPS/REcat Class II/III Transposable Elements',
                    'MIPS/REcat Other'
                ];
            },
        },
    );
    my @groups = (
        {   'type'               => 'Gene',
            'label'              => 'BAC Genes',
            'colour'             => 'navy',
            'feature_method'     => 'get_all_Genes_by_type',
            'feature_parameters' => 'genes',
        },
        {   'type'               => 'EST',
            'colour'             => 'gold2',
            'feature_parameters' => 'alignments',
        },
        {   'type'               => 'FST',
            'colour'             => 'darkorchid2',
            'feature_parameters' => 'alignments',
        },
        {   'type'               => 'GSS',
            'colour'             => 'lightcoral',
            'feature_parameters' => 'alignments',
        },
        {   'type'               => 'Array',
            'colour'             => 'seagreen4',
            'feature_parameters' => 'alignments',
        },
        {   'type'               => 'Repeat',
            'colour'             => 'firebrick',
            'feature_method'     => 'get_all_RepeatFeatures',
            'argument_index'     => 1,
            'feature_parameters' => 'repeats',
        }
    );

    for my $group (@groups) {
        my $label = $group->{label} || "$group->{type} Features";
        my $params = $dna_feature_parameters{ $group->{feature_parameters} };
        $self->add_track(
            "$group->{type}_alignments",
            'label'           => $label,
            'on'              => "on",
            'pos'             => $POS++,
            'col'             => $group->{colour},
            'lab'             => 'black',
            'width'           => 12,
            'height'          => 12,
            'glyphset'        => 'cyto_dna_alignments',
            'subtype_labeler' => $params->{subtype_labeler},
            'feature_method'  => $group->{feature_method} || undef,
            'str'             => 'r',
            'dep'             => '9999999',
            'threshold'       => '2e6',
            'sets'            => $params->{grouper}->($group->{type}),
            '_menu'           => 'features',
            '_menu_caption'   => $label,
            'argument_index'  => $group->{argument_index} || undef,
        );
    }
}

sub ADD_ALL_MARKER_FEATURES {
    my $self = shift;
    my $POS = shift || 2290;

    # Only add once!
    if ($self->{_ADDED_ALL_MARKER_FEATURES}) { return $POS }
    my %logic_names;
    for my $sp (@{ $self->{'species_defs'}->ENSEMBL_SPECIES }) {
        my $anal_adapt = $self->registry->get_adaptor($sp, 'core', 'Analysis')
            || next;
        my @analyses
            = @{ $anal_adapt->fetch_all_by_feature_class('MarkerFeature') };
        map { $logic_names{ $_->logic_name }++ } @analyses;
    }
    my $call = 'get_all_MarkerFeatures';

    my %captions
        = map { $_ => $self->logic_name_label($_) } keys %logic_names;

    # hard-coded order first, then default to alpha on captions
    my %order = (
        'core_bin_marker' => 5,
        'Maize_marker'    => 10,
        'maize-overgos'   => 20,
        'overgo-ap'       => 30,
        'overgo-dupont'   => 40,
        'RFLP'            => 50,
        'PCR'             => 60,
        'ssr_marker'      => 65,
        'CentA-int7'      => 70,
        'CentA-ltr'       => 80,
        'Telo'            => 90,
        'Chloro'          => 100,
        'Ribo'            => 110,
        'Knob'            => 120,
        'Mito'            => 130,
        'p-CentC'         => 140,
    );

    my @ordered_tracks = sort {
        ($order{$a} ||= 9999) <=> ($order{$b} ||= 9999)
            or $captions{$a} cmp $captions{$b}
    } keys %logic_names;

    my %extra = ();
    my $type  = 'Marker';

    my @analyses = grep {
               $_ !~ /marker_marker/
            && $_ ne 'core_bin_marker'
            && $_ ne 'electronic_ssr'
    } @ordered_tracks;

    $self->add_new_track_genericmatch(
        'core_bin_marker',
        'Core Markers',
        'darkgreen',
        $POS++,
        'str'     => 'r',
        'dep'     => 20,
        'labels'  => 0,
        'on'      => 'on',
        'URL_KEY' => 'GRAMENE_MARKER',
        'ZMENU'   => [
            "$type: ###ID###",
            '01:MaizeGDB'                    => '###MAIZEGDB_MARKER###',
            '02:Gramene Marker Details'      => '###HREF###',
            '03:Highlight associated clones' => '###HIGHLIGHT###',
        ],
        'CALL'      => $call,
        'THRESHOLD' => 80,
        '_menu'     => 'features',
        %extra,
    );
    $self->add_new_track_genericmatch(
        'electronic_ssr',
        'Electronic SSRs',
        'navyblue',
        $POS++,
        'str'     => 'r',
        'dep'     => 20,
        'labels'  => 0,
        'on'      => 'on',
        'URL_KEY' => 'GRAMENE_MARKER',
        'ZMENU'   => [
            "SSR: ###ID###",
            '01:MaizeGDB'                    => '###MAIZEGDB_MARKER###',
            '02:Gramene Marker Details'      => '###HREF###',
            '03:Highlight associated clones' => '###HIGHLIGHT###',
        ],
        'CALL'      => $call,
        'THRESHOLD' => 80,
        '_menu'     => 'features',
        %extra,
    );
    my $colours = +{
        'Maize_marker'  => 'seagreen2',
        'maize-overgos' => 'purple1',
        'overgo-ap'     => 'maroon2',
        'overgo-dupont' => 'coral1',
        'RFLP'          => 'darkgoldenrod1',
        'PCR'           => 'goldenrod2',
        'ssr_marker'    => 'darkolivegreen3',
        'CentA-int7'    => 'palegreen3',
        'CentA-ltr'     => 'royalblue1',
        'Telo'          => 'mediumpurple3',
        'Chloro'        => 'gold2',
        'Ribo'          => 'violetred2',
        'Knob'          => 'blueviolet',
        'Mito'          => 'blue',
        'p-CentC'       => 'blue',
    };
    $self->add_new_track_genericmatch(
        'Markers',
        'Hybridized Markers',
        'black',
        $POS++,
        'str'     => 'r',
        'dep'     => 20,
        'labels'  => 0,
        'on'      => 'on',
        'URL_KEY' => 'GRAMENE_MARKER',
        'ZMENU'   => [
            "$type: ###ID###",
            '01:MaizeGDB'                    => '###MAIZEGDB_MARKER###',
            '02:Gramene Marker Details'      => '###HREF###',
            '03:Highlight associated clones' => '###HIGHLIGHT###',
        ],
        'CALL'      => $call,
        'THRESHOLD' => 80,
        'FEATURES'  => join(' ', @analyses),
        '_menu'     => 'features',
        'colours'   => $colours,
        %extra,
    );

    $self->{_ADDED_ALL_MARKER_FEATURES}++;
    return $POS;
}

1;
