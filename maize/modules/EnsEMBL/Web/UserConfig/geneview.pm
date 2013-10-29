package EnsEMBL::Web::UserConfig::geneview;
use strict;
use EnsEMBL::Web::UserConfig;
use vars qw(@ISA);
@ISA = qw(EnsEMBL::Web::UserConfig);

sub init {
    my ($self) = @_;
    $self->{'_userdatatype_ID'} = 11;
    $self->{'fakecore'}         = 1;

    $self->{'general'}->{'geneview'} = {
        '_artefacts' =>
            [qw(ruler regulatory_search_regions regulatory_regions)],
        '_options'  => [qw(pos col known unknown)],
        'fakecore'  => 1,
        '_settings' => {
            'show_labels'                    => 'yes',
            'show_buttons'                   => 'no',
            'width'                          => 500,
            'opt_zclick'                     => 1,
            'show_empty_tracks'              => 'yes',
            'show_empty_tracks'              => 'yes',
            'bgcolor'                        => 'background1',
            'bgcolour1'                      => 'background1',
            'bgcolour2'                      => 'background1',
            'opt_protein_coding'             => 1,
            'opt_protein_coding_unsupported' => 1,
            'opt_transposon_pseudogene'      => 1,
        },
        'ruler' => {
            'on'  => 'on',
            'str' => 'r',
            'pos' => '10',
            'col' => 'black',
        },

        # col is for colours. Not needed here as overwritten in Glyphset
        'regulatory_regions' => {
            'on'        => "off",
            'pos'       => '12',
            'str'       => 'b',
            'available' => 'database_tables ENSEMBL_DB.regulatory_feature',
        },

        'regulatory_search_regions' => {
            'on'  => "off",
            'pos' => '13',
            'str' => 'b',
            'available' =>
                'database_tables ENSEMBL_DB.regulatory_search_region',
        },
    };
    $self->ADD_ALL_TRANSCRIPTS(0, 'on' => 'off');
    $self->ADD_ALL_PREDICTIONTRANSCRIPTS(0, 'on' => 'off');
}
1;
