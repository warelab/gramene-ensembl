package EnsEMBL::Web::UserConfig;

use strict;
use Data::Dumper;
use Storable qw(nfreeze freeze thaw);
use Sanger::Graphics::TextHelper;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Exception qw(throw);

use Readonly;

Readonly my @BIOTYPE_COLORS => (
    'KNOWN'         => [ 'rust',     'Known Gene' ],
    'NOVEL'         => [ 'navy',     'Novel Gene' ],
    'PSEUDOGENE'    => [ 'darkgrey', 'Pseudogene' ],
    'KNOWN_HOMOLOG' => [ 'maroon',   'Known Gene, has Ortholog' ],
    'NOVEL_HOMOLOG' => [ 'navy',     'Novel Gene, has Ortholog' ],
    'TRANSPOSON_PSEUDOGENE' =>
        [ 'maroon', 'Similarity to transposable elements' ],
    'PROTEIN_CODING' => [ 'navy', 'Similarity to known proteins' ],
    'PROTEIN_CODING_UNSUPPORTED' =>
        [ 'black', 'No similarity to known proteins' ],
    'CORRUPTED_TRANSLATION' => [ 'gray60', 'Problematic Translation' ],
);

#########
# 'general' settings contain defaults.
# 'user' settings are restored from cookie if available
# 'general' settings are overridden by 'user' settings
#

sub new {
    my $class   = shift;
    my $adaptor = shift;
    my $type    = $class =~ /([^:]+)$/ ? $1 : $class;
    my $style   = $adaptor->get_species_defs->ENSEMBL_STYLE || {};
    my $self    = {
        '_colourmap' => $adaptor->colourmap,
        '_font_face' => $style->{GRAPHIC_FONT} || 'Arial',
        '_font_size' => ($style->{GRAPHIC_FONTSIZE} * $style->{GRAPHIC_LABEL})
            || 20,
        '_texthelper'  => new Sanger::Graphics::TextHelper,
        '_db'          => $adaptor->get_adaptor,
        'type'         => $type,
        'species'      => $ENV{'ENSEMBL_SPECIES'} || '',
        'species_defs' => $adaptor->get_species_defs,
        'exturl'       => $adaptor->exturl,
        'general'      => {},
        'user'         => {},
        '_managers'        => {},       # contains list of added features....
        '_useradded'       => {},       # contains list of added features....
        '_userdatatype_ID' => 0,
        '_r'               => undef,    # $adaptor->{'r'} || undef,
        'no_load'          => undef,
        'storable'         => 1,
        'altered'          => 0
    };

    bless($self, $class);

    ########## init sets up defaults in $self->{'general'}
    $self->init()          if ($self->can('init'));
    $self->das_sources(@_) if (@_);                   # we have das sources!!

    ########## load sets up user prefs in $self->{'user'}
    #  $self->load() unless(defined $self->{'no_load'});
    return $self;
}

sub storable : lvalue {
### a
### Set whether this ScriptConfig is changeable by the User, and hence needs to
### access the database to set storable do $script_config->storable = 1; in SC code...
    $_[0]->{'storable'};
}

sub altered : lvalue {
### a
### Set to one if the configuration has been updated...
    $_[0]->{'altered'};
}

sub registry {
    my $self = shift;
    $self->{_registry} ||= "Bio::EnsEMBL::Registry";
}

sub TRIM {
    return sub { return $_[0] =~ /(^[^\.]+)\./ ? $1 : $_[0] };
}

sub update_config_from_parameter {
    my ($self, $string) = @_;
    my @array = split /\|/, $string;
    shift @array;
    return unless @array;
    foreach (@array) {
        my ($key, $value) = /^(.*):(.*)$/;
        if ($key =~ /bump_(.*)/) {
            $self->set($1, 'compact', $value eq 'on' ? 0 : 1);
        } elsif ($key eq 'imagemap' || $key =~ /^opt_/) {
            $self->set('_settings', $key, $value eq 'on' ? 1 : 0);
        } elsif ($key =~ /managed_(.*)/) {
            $self->set($key, 'on', $value, 1);
        } else {
            $self->set($key, 'on', $value);
        }
    }
    $self->save();
}

sub add_track {
    my ($self, $code, %pars) = @_;
    ## Create drop down menu entry
    my $type_config = $self->{'general'}->{ $self->{'type'} };
    if ($pars{'_menu'}) {
        $type_config->{'_settings'}{ $pars{'_menu'} } ||= [];
        push(
            @{ $type_config->{'_settings'}{ $pars{'_menu'} } },
            [ $code, $pars{'_menu_caption'} || $pars{'caption'} ]
        );
        delete $pars{'_menu'};
        delete $pars{'_menu_caption'};
    }
    push @{ $type_config->{'_artefacts'} }, $code;
    $type_config->{$code} = {%pars};
    ## Create configuration entry....
}

sub add_GSV_protein_domain_track {
    my ($self, $code, $text_label, $pos, %pars) = @_;
    $self->add_track(
        $code,
        'on'         => 'on',
        'pos'        => $pos,
        'glyphset'   => 'GSV_generic_domain',
        '_menu'      => 'features',
        'available'  => "features $code",
        'logic_name' => $code,
        'caption'    => $text_label,
        'dep'        => 20,
        'url_key'    => uc($code),
        'colours' => { $self->{'_colourmap'}->colourSet('protein_features') },
        %pars
    );
}

sub add_protein_domain_track {
    my ($self, $code, $text_label, $pos, %pars) = @_;
    $self->add_track(
        $code,
        'on'         => 'on',
        'pos'        => $pos,
        'glyphset'   => 'P_domain',
        '_menu'      => 'features',
        'available'  => "features $code",
        'logic_name' => $code,
        'caption'    => $text_label,
        'dep'        => 20,
        'url_key'    => uc($code),
        'colours' => { $self->{'_colourmap'}->colourSet('protein_features') },
        %pars
    );
}

sub add_protein_feature_track {
    my ($self, $code, $text_label, $pos, %pars) = @_;
    $self->add_track(
        $code,
        'on'         => 'on',
        'pos'        => $pos,
        'glyphset'   => 'P_feature',
        '_menu'      => 'features',
        'available'  => "features $code",
        'logic_name' => $code,
        'caption'    => $text_label,
        'colours' => { $self->{'_colourmap'}->colourSet('protein_features') },
        %pars
    );
}

sub add_new_simple_track {
    my ($self, $code, $text_label, $colour, $pos, %pars) = @_;
    $self->add_track(
        $code,
        'on'        => 'off',
        'pos'       => $pos,
        'col'       => $colour,
        'glyphset'  => 'generic_simplest',
        '_menu'     => 'features',
        'available' => "features $code",
        'str'       => 'r',
        'label'     => $text_label,
        'caption'   => $text_label,
        'code'      => $code,
        %pars
    );
}

sub add_new_synteny_track {
    my ($self, $species, $short, $pos) = @_;
    $self->add_track(
        "synteny_$species",
        "_menu"     => 'compara',
        'height'    => 4,
        'glyphset'  => "generic_synteny",
        'label'     => "$short synteny",
        'caption'   => "$short synteny",
        'species'   => $species,
        'available' => "multi SYNTENY|$species",
        'on'        => 'on',
        'pos'       => $pos,
        'str'       => 'f',
        'dep'       => 20
    );
}

sub add_new_track_transcript {
    my ($self, $code, $text_label, $colours, $pos, %pars) = @_;
    my $available = $pars{'available'} || "features $code";
    delete($pars{'available'});
    $self->add_track(
        $code . "_transcript",
        '_menu'         => 'features',
        'on'            => 'on',
        'colours'       => { $self->{'_colourmap'}->colourSet($colours) },
        'colour_set'    => $colours,
        'pos'           => $pos,
        'str'           => 'b',
        'db'            => 'core',
        'logic_name'    => $code,
        'compact'       => 0,
        'join'          => 0,
        'join_x'        => -10,
        'join_col'      => 'blue',
        'track_label'   => $text_label,
        'label'         => $text_label,
        'caption'       => $text_label,
        'available'     => $available,
        'zmenu_caption' => $text_label,
        'author'        => $pars{'author'},
        'glyphset'      => $pars{'glyph'},
        %pars
    );
}

sub add_new_track_generictranscript {
    my ($self, $code, $text_label, $colour, $pos, %pars) = @_;
    my $available = $pars{'available'} || "features $code";
    delete($pars{'available'});
    $self->add_track(
        $text_label,
        'glyphset'    => 'generic_transcript',
        'LOGIC_NAME'  => $code,
        '_menu'       => 'features',
        'on'          => 'on',
        'col'         => $colour,
        'pos'         => $pos,
        'str'         => 'b',
        'hi'          => 'highlight1',
        'compact'     => 1,
        'track_label' => $text_label,
        'caption'     => $text_label,
        'available'   => $available,
        %pars
    );
}

sub add_new_track_predictiontranscript {
    my ($self, $code, $text_label, $colour, $pos, $additional_zmenu, %pars)
        = @_;
    $self->add_track(
        $code,
        'glyphset'         => 'prediction_transcript',
        'LOGIC_NAME'       => $code,
        '_menu'            => 'features',
        'on'               => 'on',
        'col'              => $colour,
        'pos'              => $pos,
        'str'              => 'b',
        'hi'               => 'highlight1',
        'compact'          => 0,
        'track_label'      => $text_label,
        'caption'          => $text_label,
        'available'        => "features $code",
        'ADDITIONAL_ZMENU' => $additional_zmenu || {},
        %pars
    );
}

sub add_new_track_gene {
    my ($self, $code, $text_label, $colours, $pos, %pars) = @_;
    $self->add_track(
        "gene_$code",
        '_menu'                => 'features',
        'on'                   => 'on',
        'colour_set'           => $colours,
        'gene_col'             => $code,
        'pos'                  => $pos,
        'glyphset'             => 'generic_gene',
        'threshold'            => 2e6,
        'navigation_threshold' => 1e4,
        'navigation'           => 'on',
        'label_threshold'      => 1e-100,             #1e4,
        'database'             => '',
        'logic_name'           => $code,
        'available'            => "features $code",
        'caption'              => $text_label,
        'track_label'          => $text_label,
        'label'                => $text_label,
        %pars
    );
}

sub add_new_track_genericmatch {
    my ($self, $code, $text_label, $colour, $pos, %pars) = @_;
    $self->add_track(
        $code,
        '_menu'      => 'features',
        'on'         => 'off',
        'dep'        => '6',
        'str'        => 'b',
        'compact'    => 0,
        'glyphset'   => 'generic_match',
        'col'        => $colour,
        'pos'        => $pos,
        'available'  => "features $code",
        'caption'    => $text_label,
        'TEXT_LABEL' => $text_label,
        'FEATURES'   => $code,
        %pars
    );
}

sub add_new_track_cdna {
    my ($self, $code, $text_label, $pos, %pars) = @_;
    $self->add_track(
        $code,
        '_menu'      => 'features',
        'on'         => "off",
        'colour_set' => 'cdna',
        'dep'        => '6',
        'str'        => 'b',
        'compact'    => 0,
        'glyphset'   => 'generic_match',
        'SUBTYPE'    => sub {
            $_[0] =~ /^NM_/ ? 'refseq'
                : ($_[0] =~ /(RO|ZX|PX|ZA|PL)\d{5}[A-Z]\d{2}/ ? 'riken'
                : 'default');
        },
        'URL_KEY' => {
            'refseq'           => 'REFSEQ',
            'riken'            => 'RIKEN',
            'default'          => 'EMBL',
            'genoscope_ecotig' => 'TETRAODON_ECOTIG',
            'genoscope'        => 'TETRAODON_CDM'
        },
        'ID'    => { 'refseq' => TRIM },
        'ZMENU' => {
            'refseq' => [ '###ID###', "REFSEQ: ###ID###" => '###HREF###' ],
            'riken'  => [ '###ID###', "RIKEN:  ###ID###" => '###HREF###' ],
            'genoscope_ecotig' =>
                [ '###ID###', "Genoscope Ecotig:  ###ID###" => '###HREF###' ],
            'genoscope' =>
                [ '###ID###', "Genoscope:  ###ID###" => '###HREF###' ],
            'default' => [ '###ID###', "EMBL:   ###ID###" => '###HREF###' ],
        },
        'pos'        => $pos,
        'available'  => "features $code",
        'caption'    => $text_label,
        'TEXT_LABEL' => $text_label,
        'FEATURES'   => $code,
        %pars
    );
}

sub add_new_track_mrna {
    my ($self, $code, $text_label, $pos, %pars) = @_;
    $self->add_track(
        $code,
        '_menu'      => 'features',
        'on'         => "off",
        'colour_set' => 'mrna',
        'dep'        => '6',
        'str'        => 'b',
        'compact'    => 0,
        'glyphset'   => 'generic_match',
        'SUBTYPE'    => sub { return 'default'; },
        'URL_KEY'    => 'EMBL',
        'ZMENU'      => [ '###ID###', "EMBL: ###ID###" => '###HREF###' ],
        'pos'        => $pos,
        'available'  => "features $code",
        'caption'    => $text_label,
        'TEXT_LABEL' => $text_label,
        'FEATURES'   => $code,
        %pars
    );
}

sub add_new_track_rna {
    my ($self, $code, $text_label, $pos, %pars) = @_;
    $self->add_track(
        $code,
        '_menu'      => 'features',
        'on'         => "off",
        'colour_set' => 'rna',
        'dep'        => '6',
        'str'        => 'b',
        'compact'    => 0,
        'glyphset'   => 'generic_match',
        'SUBTYPE'    => $code,
        'URL_KEY'    => uc($code),
        'ZMENU' => [ '###ID###', "$text_label: ###ID###" => '###HREF###' ],
        'pos'   => $pos,
        'available'  => "features $code",
        'caption'    => $text_label,
        'TEXT_LABEL' => $text_label,
        'FEATURES'   => $code,
        %pars
    );
}

sub add_new_track_protein {
    my ($self, $code, $text_label, $pos, %pars) = @_;
    $self->add_track(
        $code,
        '_menu'      => 'features',
        'on'         => "off",
        'colour_set' => 'protein',
        'CALL'       => 'get_all_ProteinAlignFeatures',
        'dep'        => '6',
        'str'        => 'b',
        'compact'    => 0,
        'glyphset'   => 'generic_match',
        'SUBTYPE'    => sub { $_[0] =~ /^NP_/ ? 'refseq' : 'default' },
        'URL_KEY'    => 'SRS_PROTEIN',
        'ID'         => TRIM,
        'LABEL'      => TRIM,
        'ZMENU' => [ '###ID###', 'Protein homology ###ID###', '###HREF###' ],
        'pos'   => $pos,
        'available'  => "features $code",
        'caption'    => $text_label,
        'TEXT_LABEL' => $text_label,
        'FEATURES'   => $code,
        %pars
    );
}

sub add_new_track_est {
    my ($self, $code, $text_label, $pos, %pars) = @_;
    $self->add_track(
        $code,
        '_menu'      => 'features',
        'on'         => "off",
        'colour_set' => 'est',
        'dep'        => '6',
        'str'        => 'b',
        'compact'    => 0,
        'glyphset'   => 'generic_match',
        'SUBTYPE'    => sub { $_[0] =~ /^BX/ ? 'genoscope' : 'default' },
        'URL_KEY'    => 'EMBL',
        'pos'        => $pos,
        'available'  => "features $code",
        'caption'    => $text_label,
        'TEXT_LABEL' => $text_label,
        'FEATURES'   => $code,
        'ZMENU' => [ 'EST', "EST: ###ID###" => '###HREF###' ],
        %pars
    );
}

sub add_new_track_est_protein {
    my ($self, $code, $text_label, $pos, %pars) = @_;
    $self->add_track(
        $code,
        '_menu'      => 'features',
        'on'         => "off",
        'colour_set' => 'est',
        'dep'        => '6',
        'str'        => 'b',
        'compact'    => 0,
        'glyphset'   => 'generic_match',
        'CALL'       => 'get_all_ProteinAlignFeatures',
        'SUBTYPE'    => sub { $_[0] =~ /^BX/ ? 'genoscope' : 'default' },
        'URL_KEY'    => 'EMBL',
        'pos'        => $pos,
        'available'  => "features $code",
        'caption'    => $text_label,
        'TEXT_LABEL' => $text_label,
        'FEATURES'   => $code,
        'ZMENU' => [ 'EST', "EST: ###ID###" => '###HREF###' ],
        %pars
    );
}

sub add_clone_track {
    my ($self, $code, $track_label, $pos, %pars) = @_;
    $self->add_track(
        "cloneset_$code",
        '_menu'             => 'options',
        'on'                => 'off',
        'dep'               => 9999,
        'str'               => 'r',
        'glyphset'          => 'generic_clone',
        'pos'               => $pos,
        'navigation'        => 'on',
        'outline_threshold' => '350000',
        'colour_set'        => 'clones',
        'FEATURES'          => $code,
        'label'             => $track_label,
        'caption'           => $track_label,
        'available'         => 'features MAPSET_' . uc($code),
        'threshold_array'   => {
            100000 => { 'navigation' => 'off', 'height' => 4 },
            %{ $pars{'thresholds'} || {} }
        },
        %pars,
    );
}

sub add_new_track_Pgeneric_match {

    # Adds a new track to the protein feature (ProtView) image
    my ($self, $logic_name, $track_label, $colour, $pos, $compact, %pars)
        = @_;

    $self->add_track(
        $logic_name,
        'glyphset'    => 'Pgeneric_match',
        'available'   => 1,                    #"features $logic_name",
        'on'          => "on",
        'pos'         => $pos,
        'dep'         => '6',
        'col'         => $colour || 'black',
        'track_label' => $track_label,
        'caption'     => $track_label,
        'LOGIC_NAME'  => $logic_name,
        'URL_KEY'     => uc($logic_name),
        'compact'     => $compact,
        %pars
    );
}

sub set_species {
    my $self = shift;
    $self->{'species'} = shift;
}

sub get_user_settings {
    my $self = shift;
    return $self->{'user'};
}

sub artefacts {
    my $self = shift;
    return @{ $self->{'general'}->{ $self->{'type'} }->{'_artefacts'} || [] };
}

sub remove_artefacts {
    my $self = shift;
    my %artefacts = map { ($_, 1) } @_;
    @{ $self->{'general'}->{ $self->{'type'} }->{'_artefacts'} }
        = grep { !$artefacts{$_} } $self->subsections();
}

sub add_artefacts {
    my $self = shift;
    $self->_set($_, 'on', 'on') foreach @_;
    push @{ $self->{'general'}->{ $self->{'type'} }->{'_artefacts'} }, @_;
}

# add general and artefact settings
sub add_settings {
    my $self     = shift;
    my $settings = shift;
    foreach (keys %{$settings}) {
        $self->{'general'}->{ $self->{'type'} }->{$_} = $settings->{$_};
    }
}

sub turn_on {
    my $self = shift;
    $self->_set($_, 'on', 'on') foreach (@_ ? @_ : $self->subsections(1));
}

sub turn_off {
    my $self = shift;
    $self->_set($_, 'on', 'off') foreach (@_ ? @_ : $self->subsections(1));
}

sub _set {
    my ($self, $entry, $key, $value) = @_;
    $self->{'general'}->{ $self->{'type'} }->{$entry}->{$key} = $value;
}

sub load {
    my ($self) = @_;
    my @caller = caller();
    print STDERR
     "[$caller[0] Line $caller[2]] UserConfig->load - Deprecated call now handled by session\n";
    return;
    if ($self->{'_db'}) {
        my $TEMP = $self->{'_db'}
            ->getConfigByName($ENV{'ENSEMBL_FIRSTSESSION'}, $self->{'type'});
        eval { $self->{'user'} = Storable::thaw($TEMP) if $TEMP; };
    }
    return;
}

sub save {
    my ($self) = @_;
    my @caller = caller();
    print STDERR
     "[$caller[0] Line $caller[2]] UserConfig->save - Deprecated call now handled by session\n";
    return;
    $self->{'_db'}->setConfigByName(
        $self->{'_r'},   $ENV{'ENSEMBL_FIRSTSESSION'},
        $self->{'type'}, &Storable::nfreeze($self->{'user'})
    ) if ($self->{'_db'});
    return;
}

sub reset {
    my ($self) = @_;
    my $script = $self->script();
    $self->{'user'}->{$script} = {};
    $self->altered = 1;
    return;
}

sub reset_subsection {
    my ($self, $subsection) = @_;
    my $script = $self->script();
    return unless (defined $subsection);

    $self->{'user'}->{$script}->{$subsection} = {};
    $self->altered = 1;
    return;
}

sub dump {
    my ($self) = @_;
    print STDERR Dumper($self);
}

sub script {
    my ($self) = @_;
    my @keys = keys %{ $self->{'general'} };
    return $keys[0];
}

#########
# return artefacts on scripts
#
sub subsections {
    my ($self, $flag) = @_;
    my @keys;
    @keys = grep {/^managed_/} keys %{ $self->{'user'} } if $flag == 1;
    return @{ $self->{'general'}->{ $self->script }->{'_artefacts'} }, @keys;
}

#########
# return available artefacts on scripts
#
sub get_available_artefacts {
    my ($self) = @_;

    # Loop for all artefacts
    my @available_artefacts;
    foreach ($self->subsections()) {

        # Test availability
        push(@available_artefacts, $_) if $self->is_available_artefact($_);
    }

    # Return available only
    return (@available_artefacts);
}

sub required_databases {
    my $self = shift;
    my %databases;
    foreach my $a ($self->get_available_artefacts) {
        next unless $self->get($a, 'on') eq 'on';

        # get databases based on 'available' condition
        my @test = split(' ', $self->get($a, 'available'));
        if ($test[0] eq 'database_tables') {
            my ($database, $table) = split('\.', $test[1]);
            $database = lc($1) if $database =~ /^ENSEMBL_(.*)$/;
            $databases{$database} = 1;
        } elsif ($test[0] eq 'database_features') {
            my ($database, $logic_name) = split /\./, $test[1];
            $database = lc($1) if $database =~ /^ENSEMBL_(.*)$/;
            $databases{$database} = 1;
        } elsif ($test[0] eq 'databases') {
            my $database = $test[1];
            $database = lc($1) if $database =~ /^ENSEMBL_(.*)$/;
            $databases{$database} = 1;
        } elsif ($test[0] eq 'multi') {
            $databases{'compara'} = 1;
        }

        # get additional configured databases
        map { $databases{$_} = 1 } split(/,/, $self->get($a, 'databases'));
    }
    return keys %databases;
}

########
# tests whether a given artifact is available.
# data availability test for a feature is defined in the
# appropriate WebUserConfig file
# IN: self, artifact
# OUT: 999 (no test found)
#      1   (data available)
#      0   (test failed)
sub is_available_artefact {
    my $self     = shift;
    my $artefact = shift || return undef();
    my $DEBUG    = shift;
    my $settings = $self->values($artefact);
    return 0 unless $settings && %$settings;
    return $self->_is_available_artefact($settings->{available});
}

sub _is_available_artefact {
    my $self = shift;
    return $self->{'species_defs'}
        ->_is_available_artefact($self->{'species'}, @_);
}

#########
# return a list of the available options for this set of artefacts
#
sub options {
    my ($self) = @_;
    my $script = $self->script();
    return @{ $self->{'general'}->{$script}->{'_options'} };
}

sub is_setting {
    my ($self, $key) = @_;
    my $script = $self->script();
    return exists $self->{'general'}{$script}{'_settings'}{$key};
}
#########
# return a hashref of settings (user XOR general) for artefacts on scripts
#
sub values {
    my ($self, $subsection) = @_;
    my $userref;
    my $genref;
    my $hashref;

    my $script = $self->script();
    return {} unless (defined $self->{'general'}->{$script});
    return {} unless (defined $self->{'general'}->{$script}->{$subsection});

    $userref = $self->{'user'}->{$script}->{$subsection};
    $genref  = $self->{'general'}->{$script}->{$subsection};

    for my $key (keys %{$genref}) {
        $$hashref{$key} = $$userref{$key} || $$genref{$key};
    }
    return $hashref;
}

sub canset {
    my ($self, $subsection, $key) = @_;
    my $script = $self->script();

    return 1 if ($self->useraddedsource($subsection));
    return 1
        if (defined $self->{'general'}->{$script}->{$subsection}->{$key});
    return undef;
}

sub useraddedsource {
    my ($self, $subsection) = @_;
    my $useradded = $self->{'_useradded'};
    return exists $useradded->{$subsection};
}

sub set {
    my ($self, $subsection, $key, $value, $force) = @_;
    my $script = $self->script();
    return unless (defined $key && defined $script && defined $subsection);
    if ($force == 1) {
        $self->{'user'}->{$script}->{$subsection} ||= {};
    } else {
        return unless (defined $self->{'general'}->{$script});
        return unless (defined $self->{'general'}->{$script}->{$subsection});
        return
            unless (
            defined $self->{'general'}->{$script}->{$subsection}->{$key});
    }
    my ($package, $filename, $line) = caller;
    return if $self->{'user'}->{$script}->{$subsection}->{$key} eq $value;
    $self->altered = 1;
    $self->{'user'}->{$script}->{$subsection}->{$key} = $value;
}

sub get {
    my ($self, $subsection, $key) = @_;
    my $script = $self->script();

    return unless (defined $key && defined $script && defined $subsection);
    my $user_pref = undef;
    if (   defined $self->{'user'}->{$script}
        && defined $self->{'user'}->{$script}->{$subsection})
    {
        $user_pref = $self->{'user'}->{$script}->{$subsection}->{$key};
    }
    return $user_pref if (defined $user_pref);

    return unless (defined $self->{'general'}->{$script});
    return unless (defined $self->{'general'}->{$script}->{$subsection});
    return
        unless (defined $self->{'general'}->{$script}->{$subsection}->{$key});

    my $default = $self->{'general'}->{$script}->{$subsection}->{$key};
    return $default;
}

sub species_defs { return $_[0]->{'species_defs'}; }

sub colourmap {
    my ($self) = @_;
    return $self->{'_colourmap'};
}

sub image_width {
    my ($self, $script) = @_;
    return $self->{'panel_width'} || $self->get('_settings', 'width');
}

sub image_height {
    my ($self, $height) = @_;
    $self->{'_height'} = $height if (defined $height);
    return $self->{'_height'};
}

sub bgcolor {
    my ($self, $script) = @_;
    return $self->get('_settings', 'bgcolor');
}

sub bgcolour {
    my ($self, $script) = @_;
    return $self->bgcolor($script);
}

sub texthelper {
    my ($self) = @_;
    return $self->{'_texthelper'};
}

sub scalex {
    my ($self, $val) = @_;
    if (defined $val) {
        $self->{'_scalex'} = $val;
        $self->{'_texthelper'}->scalex($val);
    }
    return $self->{'_scalex'};
}

sub set_width {
    my ($self, $val) = @_;
    $self->set('_settings', 'width', $val);
}

sub container_width {
    my ($self, $val) = @_;
    if (defined $val) {
        $self->{'_containerlength'} = $val;

        my $width = $self->image_width();
        $self->scalex($width / $val) if $val;
    }
    return $self->{'_containerlength'};
}

sub transform {
    my ($self) = @_;
    return $self->{'transform'};
}

#----------------------------------------------------------------------

=head2 das_sources

  Arg [1]   : Hashref - representation of das sources;
              {$key=>{label=>$label, on=>$on_or_off, ...}}
  Function  : Adds a track to the config for each das source representation
              - Adds track to '_artefacts'
              - Adds manager of type 'das'
  Returntype: Boolean
  Exceptions:
  Caller    :
  Example   :

=cut

sub das_sources {
    my ($self, $das_sources) = @_;

    $self->{'_das_offset'} ||= 2000;
    (my $NAME = ref($self)) =~ s/.*:://;
    my $cmap = $self->{'_colourmap'};

    foreach (
        sort {
            $das_sources->{$a}->{'label'} cmp $das_sources->{$b}->{'label'}
        } keys %$das_sources
        )
    {

        my $das_source = $das_sources->{$_};

        my $on = 'off';
        if ($das_source->{'on'} eq 'on') {
            $on = 'on';
        } elsif (ref($das_source->{'on'}) eq 'ARRAY') {
            foreach my $S (@{ $das_source->{'on'} }) {
                $on = 'on' if $S eq $NAME;
            }
        }

        my $col = $das_source->{'col'};
        $col = $cmap->add_hex($col) unless $cmap->is_defined($col);
        my $manager = $das_source->{'manager'} || 'das';

        $self->{'general'}->{$NAME}->{$_} = {
            'on'         => $on,
            'pos'        => $self->{'_das_offset'},
            'col'        => $col,
            'manager'    => $manager,
            'group'      => $das_source->{'group'} || 0,
            'dep'        => $das_source->{'depth'} || 0,
            'stylesheet' => $das_source->{'stylesheet'} || 'N',
            'str'        => $das_source->{'strand'} || 'b',
            'labelflag'  => $das_source->{'labelflag'} || 'N',
            'fasta'      => $das_source->{'fasta'} || [],
        };

        push @{ $self->{'general'}->{$NAME}->{'_artefacts'} }, $_;
        $self->{'_managers'}->{$manager} ||= [];
        $self->{'_das_offset'}++;
    }

    return 1;
}

sub ADD_ALL_DNA_FEATURES {

    my $self  = shift;
    my $POS   = shift || 2300;
    my %extra = @_;

    if ($self->{_ADDED_ALL_DNA_FEATURES}) { return $POS }

    my %def = map { lc($_) => 1 }    # Trackts turned on by default
        qw(maize_est maize_mrna maize_markers);

    my %colours = (
        'SWISSPROT'   => 'green',
        'RICE'        => 'limegreen',
        'MAIZE'       => 'red',
        'WHEAT'       => 'blue',
        'SUGARCANE'   => 'darkgreen',
        'RYEGRASS'    => 'darkorange',
        'SORGHUM'     => 'orange',
        'BARLEY'      => 'purple',
        'BARLEY1'     => 'purple',
        'MILLET'      => 'brown',
        'ARABIDOPSIS' => 'green',
        'BRASSICA'    => 'purple'
    );

    my %menu_types;

    $menu_types{gss} = join "|", reverse sort qw(
        HI_COT_BENNETZEN
        METH_FILT_CSHL_MCCOMBIE
        METH_FILT_TIGR
        GSS-READ_KLEIN
        MAGI_ISU
        BACEND
        METH_FILT_HI_COT_CLUSTER
        HI_COT_TIGR
        ASSEMBLY
        SEQUENCE
        ORION
    );

    $menu_types{est} = join "|", reverse sort qw(
        CDS
        MRNA
        EST
        ESTCLUSTER3P_LGBPRATT
        ESTCLUSTER_PLANTGDB
        GI
        IND_EST
        IND_CLUSTER
        CDNA_KOME
        CONSENSUS
    );

    $menu_types{array} = join "|", reverse sort qw(
        ARRAYTARGET_NSF58K
        ARRAYCONSENSUS
        MPSSTAG
        ARRAYOLIGO_NSF
        SAGETAG_MGOS
        ARRAYOLIGO
    );

    $menu_types{fst} = join "|", reverse sort qw(
        MU_INSERT
        DS_INSERT
        T_DNA_INSERT
        TRANSPOSON_INSERT_SITE
        TOS17_INSERT
        FSTTRANSPOSON
    );

    # analysis to skip;
    #my $logic_name_skip = join "|", reserver sort qw(marker_marker$ );

    # get analysis to process
    my %logic_names;
    foreach my $sp (@{ $self->{'species_defs'}->ENSEMBL_SPECIES }) {
        my $anal_adapt = $self->registry->get_adaptor($sp, 'core', 'Analysis')
            || next;
        my @analyses;
        push(
            @analyses,
            @{  $anal_adapt->fetch_all_by_feature_class('DnaAlignFeature')
                },
            @{ $anal_adapt->fetch_all_by_feature_class('MarkerFeature') },

        );
        map { $logic_names{ $_->logic_name } = $_ } @analyses;
    }

    # Note on the sort -
    # Push Rice and Swissprot to the top,
    # Then alphabetic sort by the caption (if available) or logic name.

    foreach my $logic (
        map { $_->[1] }
        sort { ($a->[0] <=> $b->[0]) || ($a->[2] cmp $b->[2]) }
        map {
            my $i = 3;
            my $l = $self->logic_name_label($_);
            $l =~ s/_/ /g;
            if    (/^swissp/i) { $i = 1 }
            elsif (/^or/i)     { $i = 2 }
            elsif (/^rice/i)   { $i = 2 }
            [ $i, $_, $l ]
        } keys %logic_names
        )
    {

        my ($species, $type) = split('_', uc($logic), 2);

        my $colour = $colours{$species} || 'black';
        my $caption = $self->logic_name_label($logic)
            || next;    # Skip if empty caption

        my @ZMENU;
        my $URL_KEY = '';
        my $ID      = '';
        my $ID1     = '';
        my $ID2     = '';
        my %EXTRA;
        my $FEATUREVIEW = 1;

    TRACK: {

            #MARKERS
            if ($type =~ /^(RFLP_MARKER|MARKER|SSR)/) {
                $EXTRA{_menu} = "markers";

                @ZMENU = ("Marker: ###ID###");

                if ($species eq 'MAIZE') {
                    $URL_KEY = 'MAIZEGDB_MARKER';
                    push(@ZMENU, '01:MaizeGDB' => '###MAIZEGDB_MARKER###');
                    push(@ZMENU, '02:CMap'     => '###CMAP_FEATURE###');
                    push @ZMENU, ("03:GenBank" => "###ENTREZ_NUCLEOTIDE###");

                } else {    #for othere species, mostly rice

                    if ($type =~ /^(RFLP_MARKER|SSR)$/) {
                        push @ZMENU, ("01:Gramene Marker" => "###MARKER###");
                        push @ZMENU,
                            ("02:View all hits" => "###MARKERVIEW###");
                        $EXTRA{'THRESHOLD'} = 0;
                        $EXTRA{'CALL'}      = 'get_all_MarkerFeatures';

                        if ($type eq 'SSR') {
                            $EXTRA{'NO_ALIGNMENT'} = 1;
                            $FEATUREVIEW = 0;
                        }
                    } else {

                        #$URL_KEY = 'ENTREZ_NUCLEOTIDE';
                        push @ZMENU,
                            (
                            "01:GenBank" => "###ENTREZ_NUCLEOTIDE###",
                            "02:Marker", "###MARKER###"
                            );
                    }
                }

                last TRACK;
            }

            my $menu_container;
            for my $menu (sort keys %menu_types) {
                if ($type =~ /$menu_types{$menu}/) {
                    $EXTRA{_menu} = "${menu}_features";
                    $menu_container = $menu;
                    last;
                }
            }

            if ($type =~ /QTL/) {
                push(@ZMENU, "01:QTL_db" => "###QTL###");
                last TRACK;
            }

            #default links
            @ZMENU = (
                "ID: ###ID###",
                "01:GenBank" => "###ENTREZ_NUCLEOTIDE###",
                "02:Marker", "###MARKER###"
            );
            $ID = sub { (split('\.', $_[0]))[0] };    # Strip version

            if ($type =~ /(METH_FILT_HI_COT_CLUSTER|HI_COT_TIGR)/) {
                push(@ZMENU, "03:TIGR" => "###TIGR_AZM_MAIZE###");
                last TRACK;
            }

            if ($type =~ /(ESTCLUSTER_PLANTGDB|MU_INSERT)/) {
                push(@ZMENU, "03:PLANTGDB" => "###PLANTGDB###");
                last TRACK;
            }

            if ($type =~ /ARRAYCONSENSUS/) {
                push(@ZMENU, "03:PLEXdb" => "###${species}PLEX###");
                $ID
                    = sub { (split(':', $_[0]))[-1] }; #strip initial consensus:SPECIES:
                last TRACK;
            }

            if ($type eq 'MPSSTAG') {
                push(@ZMENU, "03:RiceMPSS" => "###RICEMPSS###");
                $EXTRA{'NO_ALIGNMENT'} = 1;
                last TRACK;
            }

            if ($type =~ /ARRAYOLIGO_NSF/) {
                push(@ZMENU, "03:TIGR RiceArray" => "###RICEARRAY###");
                last TRACK;
            }

            if ($type =~ /^SAGETAG_MGOS/) {
                push(@ZMENU, "03:MGOS RiceSage" => "###RICESAGE###");
                last TRACK;
            }

            if ($type eq 'GI') {
                push(@ZMENU, "03:TIGR GI" => "###TIGR_GI_$species###");
                last TRACK;

            }

            if (   $type eq 'IND_EST'
                or $type eq 'IND_CLUSTER')
            {
                push(@ZMENU, "03:BGI" => "###BGI_HOME###");
                last TRACK;
            }

            if ($type eq 'CDNA_KOME') {
                push(@ZMENU, "03:KOME cDNA" => "###KOME_CDNA###");
                last TRACK;
            }

            if ($type eq 'BAC') {
                $EXTRA{labels} = 'on';
                last TRACK;
            }

            if ($type eq 'TOS17_INSERT') {
                $EXTRA{'THRESHOLD'} = 0;
                $FEATUREVIEW = 0;
                last TRACK;
            }

            if ($type eq 'CORNSENSUS') {
                push(@ZMENU, "03:MaizeGDB" => "###MAIZEGDB_OVERGO###");
                last TRACK;

            }

            if (   $logic =~ /RYEGRASS_(ASSEMBLY|SEQUENCE)/i
                || $type =~ /ORION/)
            {
                push(@ZMENU, "03:ORION" => "###ORION_CREDIT###");
                last TRACK;
            }

            if ($logic =~ /Sorghum_ESTCluster_LGBPratt/i) {
                push(@ZMENU, '03:LGB_Pratt' => '###PRATT_ESTCLUSTER###');
                last TRACK;
            }

            if ($logic =~ /MAIZE_EST/i) {
                push(@ZMENU, '03:MaizeGDB' => '###MAIZEGDB_EST###');
                last TRACK;
            }

            if ($logic =~ /Other-poaceae_est/i) {
                push(@ZMENU, '03:GrainGenes' => '###GRAINGENES_EST###');
                last TRACK;
            }

            #the rest in feature menu

        }    #end 0f Track

        push(@ZMENU, '10:View all hits' => '###FEATUREVIEW_DNA###')
            if $FEATUREVIEW;

        $self->add_new_track_genericmatch(
            $logic, $caption, $colour, $POS++,
            'str'     => 'r',
            'dep'     => 0,
            'on'      => $def{ lc($logic) } ? 'on' : 'off',
            'SUBTYPE' => $type,
            'ZMENU'   => [@ZMENU],
            'URL_KEY' => $URL_KEY,
            $ID ? ('ID' => $ID) : (),
            %EXTRA, %extra,
        );

    }

    $self->{_ADDED_ALL_DNA_FEATURES}++;
    return $POS;
}

sub ADD_ALL_MARKER_FEATURES {
    my $self = shift;
    my $POS = shift || 2290;
    return $POS;

    # Only add once!
    if ($self->{_ADDED_ALL_MARKER_FEATURES}) { return $POS }
    my %logic_names;
    foreach my $sp (@{ $self->{'species_defs'}->ENSEMBL_SPECIES }) {
        my $anal_adapt = $self->registry->get_adaptor($sp, 'core', 'Analysis')
            || next;
        my @analyses
            = @{ $anal_adapt->fetch_all_by_feature_class('MarkerFeature') };
        map { $logic_names{ $_->logic_name }++ } @analyses;
    }

    my %colours = (
        'RICE_SSR'    => 'plum3',
        'RICE_MARKER' => 'magenta',
        'TOS17'       => 'plum4'
    );
    my %types = (
        'RICE_SSR'    => 'SSR',
        'RICE_MARKER' => 'RFLP',
        'TOS17'       => 'Tos17'
    );
    my %default = (
        'RICE_SSR'    => 1,
        'RICE_MARKER' => 1,
        'TOS17'       => 1
    );
    my $call = 'get_all_MarkerFeatures';

    my %captions
        = map { $_ => $self->logic_name_label($_) } keys %logic_names;

    foreach my $logic (
        sort { $captions{$a} cmp $captions{$b} }
        keys %logic_names
        )
    {
        next if $logic =~ /marker_marker$/i;
        my $caption = $captions{$logic} || next;    # Skip if empty caption
        my $colour = $colours{ uc($logic) } || 'magenta';
        my $type   = $types{ uc($logic) }   || 'Marker';
        my %extra  = ();

        my $menu = 'markers';
        if ($logic =~ /tos17/i) { $menu = 'fst_features' }

        #warn "$logic,$caption,$type,$menu";

        if ($type eq 'SSR') { $extra{'NO_ALIGNMENT'} = 1 }

        $self->add_new_track_genericmatch(
            $logic, $caption, $colour,
            $POS++,
            'str'     => 'r',
            'dep'     => 20,
            'on'      => $default{ uc($logic) } ? 'on' : 'off',
            'URL_KEY' => 'GRAMENE_MARKER',
            'ZMENU'   => [
                "$type: ###ID###",
                '02:Gramene Marker' => '###HREF###',
                '10:View all hits'  => '###MARKERVIEW###'
            ],

            #          'labels'    => 'on',
            'CALL'      => $call,
            'THRESHOLD' => 0,
            '_menu'     => $menu,
            %extra,
        );

    }
    $self->{_ADDED_ALL_MARKER_FEATURES}++;
    return $POS;
}

sub ADD_ALL_VARIATION_TRACKS {
    my $self = shift;
    my $POS = shift || 5000;

    # Only add once!
    if ($self->{_ADDED_ALL_VARIATION_TRACKS}) { return $POS }
    my %populations;
    for my $sp (@{ $self->{'species_defs'}->ENSEMBL_SPECIES }) {
        my $pop_adapt
            = $self->registry->get_adaptor($sp, 'variation', 'Population');
        $pop_adapt || next;
        for my $pop (@{ $pop_adapt->fetch_all }) {

     # We will display a track for each SNP sample/population found in the DB.
     # Please ensure that only samples with alleles are in the DB table, i.e.;
     #  select s.sample_id, s.name, count(*)
     #  from sample s, allele a
     #  where s.sample_id=a.sample_id group by s.sample_id;
            my $pop_name   = $pop->name;
            my $track_name = $pop_name;
            $track_name =~ s/\W/_/g;
            $self->add_track(
                lc($track_name),
                'glyphset'   => 'variation',
                'on'         => 'on',          #'off',
                'bump_width' => 0,
                'dep'        => 0.1,
                'pos'        => $POS,
                'str'        => 'r',
                'colours' =>
                    { $self->{'_colourmap'}->colourSet('variation') },
                'populations' => [$pop_name],
                'track_label' => "SNPs ($pop_name)",
                'caption'     => "SNPs ($pop_name)",
                '_menu'       => 'features',
                'available'   => 'databases ENSEMBL_VARIATION'
            );
            $POS++;
        }
    }
    $self->{_ADDED_ALL_VARIATION_TRACKS}++;
    return $POS;
}

sub ADD_ALL_EST_FEATURES {
    my $self  = shift;
    my $POS   = shift || 2300;
    my %extra = @_;

    # Nothing to do for Maize
    return $POS;
}

sub ADD_ALL_CLONE_TRACKS {
    my $self = shift;
    my $POS = shift || 2500;
    $self->add_clone_track('cloneset_0_5mb', '0.5Mb clones',   $POS++, @_);
    $self->add_clone_track('cloneset_1mb',   '1Mb clones',     $POS++, @_);
    $self->add_clone_track('cloneset_30k',   '30k TPA clones', $POS++, @_);
    $self->add_clone_track('cloneset_32k',   '32k clones',     $POS++, @_);
    $self->add_clone_track('acc_bac_map',    'Acc. BAC map',   $POS++, @_);
    $self->add_clone_track(
        'bac_map', 'BAC map', $POS++,
        'threshold'  => 2e8,
        'thresholds' => { 20000 => { 'FEATURES' => 'acc_bac_map' }, },
        @_
    );

    $self->add_clone_track('bacs',              'BACs',           $POS++, @_);
    $self->add_clone_track('bacs_bands',        'Band BACs',      $POS++, @_);
    $self->add_clone_track('extra_bacs',        'Extra BACs',     $POS++, @_);
    $self->add_clone_track('tilepath_cloneset', 'Mouse Tilepath', $POS++, @_);
    $self->add_clone_track('tilepath', 'Human tilepath clones', $POS++, @_);
    $self->add_clone_track(
        'fosmid_map', 'Fosmid map', $POS++,
        'colour_set' => 'fosmids',
        @_
    );
}

sub ADD_ALL_PROTEIN_FEATURES {
    my $self  = shift;
    my $POS   = shift || 2200;
    my %extra = @_;

    # Only add once!
    if ($self->{_ADDED_ALL_PROTEIN_FEATURES}) { return $POS }

    my %logic_names;
    for my $sp (@{ $self->{'species_defs'}->ENSEMBL_SPECIES }) {
        my $anal_adapt = $self->registry->get_adaptor($sp, 'core', 'Analysis')
            || next;
        my @analyses = @{ $anal_adapt->fetch_all_by_feature_class(
                'ProteinAlignFeature') };
        map { $logic_names{ $_->logic_name }++ } @analyses;
    }

    my %captions = (
        'SPARAB'  => 'Arab_Prot_SP',
        'SPPLANT' => 'Plant_Prot_SP',
        'SPOTHER' => 'NonPlant_Prot_SP',
        'TRARAB'  => 'Arab_Prot_TrEMBL',
    );
    my %colours = (
        'SWISSPROT'   => 'green',
        'RICE'        => 'green',
        'MAIZE'       => 'red',
        'WHEAT'       => 'blue',
        'SUGARCANE'   => 'darkgreen',
        'RYEGRASS'    => 'darkorange',
        'SORGHUM'     => 'orange',
        'BARLEY'      => 'purple',
        'BARLEY1'     => 'purple',
        'MILLET'      => 'brown',
        'ARABIDOPSIS' => 'green',
        'BRASSICA'    => 'purple',
        'SPARAB'      => 'green',
        'SPPLANT'     => 'purple',
        'SPOTHER'     => 'grey',
        'TRARAB'      => 'green',
    );

    for my $logic (
        sort {
            ($captions{ uc($a) } || uc($a)) cmp($captions{ uc($b) } || uc($b))
        }
        keys %logic_names
        )
    {
        my ($species, $type) = split('_', uc($logic), 2);
        my $colour  = $colours{$species}      || 'black';
        my $caption = $captions{ uc($logic) } || $logic;
        $self->add_new_track_protein(
            $logic,
            $caption,
            $POS++,
            'str'        => 'r',
            'colour_set' => '',        #unset
            'col'        => $colour,
            %extra,
        );
    }
    $self->{_ADDED_ALL_PROTEIN_FEATURES}++;
    return $POS;
}

sub ADD_ALL_PREDICTIONTRANSCRIPTS {
    my $self  = shift;
    my $POS   = shift || 2100;
    my %extra = @_;

    # Only add once!
    if ($self->{_ADDED_ALL_PREDICTIONTRANSCRIPTS}) { return $POS }

    my %logic_names;
    for my $sp (@{ $self->{'species_defs'}->ENSEMBL_SPECIES }) {
        my $anal_adapt = $self->registry->get_adaptor($sp, 'core', 'Analysis')
            || next;
        my @analyses = @{ $anal_adapt->fetch_all_by_feature_class(
                'PredictionTranscript') };
        map { $logic_names{ $_->logic_name }++ } @analyses;
    }

    for my $logic (keys %logic_names) {
        next if $logic eq 'Fgenesh';
        my $on = 'off';
        $self->add_new_track_predictiontranscript(
            $logic,
            $self->logic_name_label($logic),
            'lightseagreen',
            $POS++,
            {},
            'on' => $on,
            %extra
        );
    }

    $self->{_ADDED_ALL_TRANSCRIPTS}++;
    return $POS;

}

sub ADD_SYNTENY_TRACKS {
    my $self = shift;
    my $POS = shift || 99900;
    for (sort @{ $self->{'species_defs'}->ENSEMBL_SPECIES }) {
        $self->add_new_synteny_track(
            $_,
            $self->{'species_defs'}->other_species($_, 'SPECIES_COMMON_NAME'),
            $POS++,
            @_
        );
    }
}

sub ADD_ALL_TRANSCRIPTS {
    my $self = shift;
    my $POS = shift || 2000;

    # Only add once!
    if ($self->{_ADDED_ALL_TRANSCRIPTS}) { return $POS }

    my %logic_names;
    for my $sp (@{ $self->{'species_defs'}->ENSEMBL_SPECIES }) {

        # Need to deal with both core and OTHERFEATURES databases
        for my $db ('OTHERFEATURES', 'core') {    # Last entry is default
                # for my $db ('OTHERFEATURES') {    # Last entry is default
            my $anal_adapt
                = $self->registry->get_adaptor($sp, $db, 'Analysis')
                || next;
            my @analyses
                = @{ $anal_adapt->fetch_all_by_feature_class('Gene') };

        # my @analyses = @{ $anal_adapt->fetch_all_by_feature_class('Gene') };
            map { $logic_names{ $_->logic_name } = $db } @analyses;
        }
    }

    for my $logic (keys %logic_names) {
        my $on = 'off';

        my ($colour, $tecolour);
        my %extra = ('db_alias' => $logic_names{$logic});

        if ($logic eq 'GeneBuilder') {
            $on     = 'off';
            $colour = 'darkgreen';
        } else {
            $on     = 'on';
            $colour = undef;    #'darkred';
            $extra{colour} = undef, $extra{colours} = {@BIOTYPE_COLORS};
        }
        if ($logic_names{$logic} ne 'core') {
            $extra{available} = "databases ENSEMBL_$logic_names{$logic}";
        }

        # print STDERR "ADDING generic_transcript $logic\n";

        $self->add_new_track_generictranscript(
            $logic,
            $self->logic_name_label($logic),
            $colour,
            $POS++,
            'on'                   => $on,
            'highlight_homologues' => '1',
            %extra,
        );
    }

    $self->{_ADDED_ALL_TRANSCRIPTS}++;
    return $POS;
}

sub ADD_ALL_OLIGO_TRACKS {
    my $self = shift;
    my $POS  = shift || 4000;
    my @AFFY = qw(
        AFFY_HG_Focus
        AFFY_HG_U133_Plus_2 AFFY_HG_U133A_2  AFFY_HG_U133A  AFFY_HG_U133B
        AFFY_HG_U95Av2      AFFY_HG_U95B     AFFY_HG_U95C   AFFY_HG_U95D   AFFY_HG_U95E
        AFFY_MG_U74Av2      AFFY_MG_U74Bv2   AFFY_MG_U74Cv2
        AFFY_Mouse430_2     AFFY_Mouse430A_2
        AFFY_RG_U34A        AFFY_RG_U34B     AFFY_RG_U34C   AFFY_Rat230_2
        AFFY_Zebrafish
    );

    for my $chipset (@AFFY) {
        (my $T  = lc($chipset)) =~ s/-/_/g;
        (my $T2 = $chipset)     =~ s/affy_/AFFY /i;
        $self->add_track(
            $T,
            'on'        => 'off',
            'pos'       => $POS++,
            'str'       => 'b',
            '_menu'     => 'features',
            'caption'   => $T2,
            'dep'       => 6,
            'col'       => 'springgreen4',
            'compact'   => 0,
            'available' => "features mapset_$T",
            'glyphset'  => 'generic_microarray',
            'FEATURES'  => $chipset,
        );
    }
    return $POS;
}

sub ADD_SIMPLE_TRACKS {
    my $self = shift;
    my $POS = shift || 7500;
    $self->add_new_simple_track('abberation_junction', 'Abberation junction',
        'red', $POS++, @_);
    $self->add_new_simple_track('enhancer', 'Enhancer', 'red', $POS++, @_);
    $self->add_new_simple_track('transcription_start_site',
        'Transcription start site',
        'red', $POS++, @_);
    $self->add_new_simple_track('regulatory_region', 'Regulatory region',
        'red', $POS++, @_);
    $self->add_new_simple_track('regulatory_search_region',
        'Regulatory search region',
        'red', $POS++, @_);
    $self->add_new_simple_track('mature_peptide', 'Mature peptide',
        'red', $POS++, @_);
    $self->add_new_simple_track('insertion_site', 'Insertion site',
        'red', $POS++, @_);
    $self->add_new_simple_track('protein_binding_site',
        'Protein binding site',
        'red', $POS++, @_);
    $self->add_new_simple_track('scaffold', 'Scaffold', 'red', $POS++, @_);
    $self->add_new_simple_track('allele',   'Allele',   'red', $POS++, @_);

#  $self->add_new_simple_track( 'RNAi',                     'RNAi',              'red', $POS++, @_ );
    $self->add_new_simple_track('fosmid', 'Fosmid', 'red', $POS++, @_);
    $self->add_new_simple_track(
        'transposable_element_insertion_site',
        'Transposable element insertion site',
        'red', $POS++, @_
    );
    $self->add_new_simple_track('transposable_element',
        'Transposable element',
        'red', $POS++, @_);
    $self->add_new_simple_track('rescue_fragment', 'Rescue fragment',
        'red', $POS++, @_);
    $self->add_new_simple_track('signal_peptide', 'Signal peptide',
        'red', $POS++, @_);
}

sub ADD_GENE_TRACKS {
    my $self = shift;
    my $POS = shift || 2000;

    # Only add once!
    if ($self->{_ADDED_ALL_GENES}) { return $POS }

    my %logic_names;
    for my $sp (@{ $self->{'species_defs'}->ENSEMBL_SPECIES }) {
        my $anal_adapt = $self->registry->get_adaptor($sp, 'core', 'Analysis')
            || next;
        my @analyses = @{ $anal_adapt->fetch_all_by_feature_class('Gene') };
        map { $logic_names{ $_->logic_name }++ } @analyses;
    }

    foreach my $logic (keys %logic_names) {

        if ($logic eq 'cDNA_arabidopsis') {next}    # Hard-coded skip
        if ($logic eq 'MinT')             {next}    # Hard-coded skip
        my $on = ($logic =~ /submitter/i) ? 'off' : 'on';
        $self->add_new_track_gene(
            $logic,
            $self->logic_name_label($logic),
            'ensembl_gene',
            $POS++,
            'colour_set' => undef,                  # Do not use builtin set
            'colours'    => {@BIOTYPE_COLORS},
            'gene_label' => sub { $_[0]->external_name || 'NOVEL' },
            'gene_col' => sub { (uc $_[0]->biotype) },
            'on'       => $on
        );
    }
    return $POS;
}

sub ADD_ALL_AS_TRANSCRIPTS {
    my $self = shift;
    my $POS = shift || 2000;
    $self->add_new_track_transcript('ensembl', 'Ensembl genes',
        'ensembl_gene', $POS++, @_);
    $self->add_new_track_transcript(
        'evega', 'Vega genes', 'vega_gene', $POS++,
        'available' => 'databases ENSEMBL_VEGA',
        @_
    );

    return $POS;
}

sub ADD_AS_GENE_TRACKS {
    my $self = shift;
    my $POS = shift || 2000;
    $self->add_new_track_gene(
        'ensembl',
        'Ensembl Genes',
        'ensembl_gene',
        $POS++,
        'gene_label' => sub {
            return $_[0]->biotype eq 'bacterial_contaminant' ? 'Bac. cont.'
                : (
                $_[0]->biotype eq 'pseudogene' ? 'Pseudogene'
                : ($_[0]->external_name || 'NOVEL')
                );
        },
        'gene_col' => sub {
            return $_[0]->biotype eq 'bacterial_contaminant' ? '_BACCOM'
                : (
                $_[0]->biotype eq 'pseudogene' ? '_PSEUDO'
                : '_' . $_[0]->external_status
                );
        },
        'logic_name' => 'ensembl psuedogene',
        @_
    );
    $self->add_new_track_gene(
        'flybase',
        'Flybase Genes',
        'flybase_gene',
        $POS++,
        'gene_label' => sub {
            return $_[0]->biotype eq 'bacterial_contaminant' ? 'Bac. cont.'
                : (
                $_[0]->biotype eq 'pseudogene' ? 'Pseudogene'
                : ($_[0]->external_name || 'NOVEL')
                );
        },
        'gene_col' => sub {
            return $_[0]->biotype eq 'bacterial_contaminant' ? '_BACCOM'
                : (
                $_[0]->biotype eq 'pseudogene' ? '_PSEUDO'
                : '_' . $_[0]->external_status
                );
        },
        'logic_name' => 'flybase psuedogene',
        @_
    );
    $self->add_new_track_gene(
        'wormbase',
        'Wormbase Genes',
        'wormbase_gene',
        $POS++,
        'gene_label' => sub {
            return $_[0]->biotype eq 'bacterial_contaminant' ? 'Bac. cont.'
                : (
                $_[0]->biotype eq 'pseudogene' ? 'Pseudogene'
                : ($_[0]->external_name || 'NOVEL')
                );
        },
        'gene_col' => sub {
            return $_[0]->biotype eq 'bacterial_contaminant' ? '_BACCOM'
                : (
                $_[0]->biotype eq 'pseudogene' ? '_PSEUDO'
                : '_' . $_[0]->external_status
                );
        },
        'logic_name' => 'wormbase psuedogene',
        @_
    );
    $self->add_new_track_gene(
        'genebuilderbeeflymosandswall',
        'Bee Genes',
        'bee_gene',
        $POS++,
        'gene_label' => sub {
            return $_[0]->biotype eq 'bacterial_contaminant' ? 'Bac. cont.'
                : (
                $_[0]->biotype eq 'pseudogene' ? 'Pseudogene'
                : ($_[0]->external_name || 'NOVEL')
                );
        },
        'gene_col' => sub {
            return $_[0]->biotype eq 'bacterial_contaminant' ? '_BACCOM'
                : (
                $_[0]->biotype eq 'pseudogene' ? '_PSEUDO'
                : '_' . $_[0]->external_status
                );
        },
        @_
    );
    $self->add_new_track_gene(
        'SGD',
        'SGD Genes',
        'sgd_gene',
        $POS++,
        'gene_label' => sub {
            return $_[0]->biotype eq 'bacterial_contaminant' ? 'Bac. cont.'
                : (
                $_[0]->biotype eq 'pseudogene' ? 'Pseudogene'
                : ($_[0]->external_name || 'NOVEL')
                );
        },
        'gene_col' => sub {
            return $_[0]->biotype eq 'bacterial_contaminant' ? '_BACCOM'
                : (
                $_[0]->biotype eq 'pseudogene' ? '_PSEUDO'
                : '_' . $_[0]->external_status
                );
        },
        @_
    );
    $self->add_new_track_gene(
        'gsten',
        'Genoscope Genes',
        'genoscope_gene',
        $POS++,
        'gene_label' => sub { return $_[0]->stable_id },
        'gene_col'   => sub {
            return $_[0]->biotype eq 'Genoscope_predicted'
                ? '_GSTEN'
                : '_HOX';
        },
        'logic_name' => 'gsten hox cyt',
        @_
    );

    $self->add_new_track_gene(
        'otter', 'Vega Genes', 'vega_gene', $POS++,
        'database'  => 'vega',
        'available' => 'databases ENSEMBL_VEGA',
        'gene_col' =>
            sub { return $_[0]->biotype . '_' . $_[0]->confidence; },
        'gene_label' => sub { $_[0]->external_name || $_[0]->stable_id; },
        @_
    );

    #for genes in Vega
    $self->add_new_track_gene(
        'vega_gene', 'Vega Genes', 'vega_gene', $POS++,
        'available'  => 'features VEGA_GENES_OTTER',
        'glyphset'   => 'vega_gene',
        'logic_name' => 'otter',
        'gene_col'   => 'vega_gene',
        @_
    );
    $self->add_new_track_gene(
        'vega_corf_gene', 'CORF Genes', 'vega_gene', $POS++,
        'available'  => 'features VEGA_GENES_OTTER_CORF',
        'glyphset'   => 'vega_gene',
        'logic_name' => 'otter_corf',
        'gene_col'   => 'vega_gene',
        @_
    );

    return $POS;
}

sub ADD_ALL_PROTEIN_FEATURE_TRACKS {
    my $self = shift;
    my $POS = shift || 2000;
    $self->add_protein_domain_track('Prints', 'PRINTS', $POS++);
    $self->add_protein_domain_track('PrositePatterns', 'Prosite patterns',
        $POS++);
    $self->add_protein_domain_track('scanprosite', 'Prosite patterns',
        $POS++);
    $self->add_protein_domain_track('PrositeProfiles', 'Prosite profiles',
        $POS++);
    $self->add_protein_domain_track('pfscan', 'Prosite profiles', $POS++);

    $self->add_protein_domain_track('Pfam',        'Pfam',            $POS++);
    $self->add_protein_domain_track('TigrFam',     'TIGRFAM',         $POS++);
    $self->add_protein_domain_track('SuperFamily', 'SUPERFAMILY',     $POS++);
    $self->add_protein_domain_track('Smart',       'SMART',           $POS++);
    $self->add_protein_domain_track('PIRS',        'PIR SuperFamily', $POS++);

    $self->add_protein_feature_track('ncoils',  'Coiled coils',     $POS++);
    $self->add_protein_feature_track('SignalP', 'Sig.Pep cleavage', $POS++);
    $self->add_protein_feature_track('Seg',     'Low complex seq',  $POS++);
    $self->add_protein_feature_track('tmhmm',   'Transmem helices', $POS++);
}

sub ADD_ALL_PROTEIN_FEATURE_TRACKS_GSV {
    my $self = shift;
    my $POS = shift || 2000;
    $self->add_GSV_protein_domain_track('Prints', 'PRINTS', $POS++);
    $self->add_GSV_protein_domain_track('PrositePatterns', 'Prosite patterns',
        $POS++);
    $self->add_GSV_protein_domain_track('scanprosite', 'Prosite patterns',
        $POS++);
    $self->add_GSV_protein_domain_track('PrositeProfiles', 'Prosite profiles',
        $POS++);
    $self->add_GSV_protein_domain_track('pfscan', 'Prosite profiles', $POS++);

    $self->add_GSV_protein_domain_track('Pfam',        'PFam',        $POS++);
    $self->add_GSV_protein_domain_track('Tigrfam',     'TIGRFAM',     $POS++);
    $self->add_GSV_protein_domain_track('Superfamily', 'SUPERFAMILY', $POS++);
    $self->add_GSV_protein_domain_track('Smart',       'SMART',       $POS++);
    $self->add_GSV_protein_domain_track('PIRSF', 'PIR SuperFamily', $POS++);
}

sub logic_name_label {
    my $self       = shift;
    my $logic_name = shift;

    my $sp = 'Zea_mays2';
    my $analysis_adaptor
        = $self->registry->get_adaptor($sp, 'core', 'Analysis');
    my $analysis_obj = $analysis_adaptor->fetch_by_logic_name($logic_name)
        if $analysis_adaptor;
    my $display_label = $analysis_obj->display_label if $analysis_obj;
    return $display_label if $display_label;

    # Can be used to disable tracks - set to empty string
    my %labels = (

        # Feature Tracks
        'SWISSPROT_TREMBL_PROTEINS' => 'Rice Protein SpTrEMBL',
        'TIGR_GENE'                 => 'Rice GeneModel TIGR',
        'TIGR'                      => 'Arabid GeneModel TIGR',
        'GENEMODEL_TIGR'            => 'Rice GeneModel TIGR',
        'SUBMITTERGENEANNOTATION'   => 'Rice GeneModel Submitted',
        'FGENESH'                   => 'Fgenesh Models',

        # EST tracks
        'RICE_MRNA'                  => 'Rice mRNA',
        'RICE_EST'                   => 'Rice EST',
        'RICE_GI'                    => 'Rice TGI ESTCl.',
        'RICE_TUG'                   => 'Rice TUG ESTCl.',
        'RICE_ESTCLUSTER_PLANTGDB'   => 'Rice PlantGDB ESTCl.',
        'RICE_IND_CLUSTER'           => 'O.in. BGI ESTCl.',
        'RICE_IND_EST'               => 'O.in. BGI ESTs',
        'RICE_JAP_CDNA_KOME'         => 'O.ja. cDNA KOME',
        'BARLEY_EST'                 => 'Barley ESTs',
        'BARLEY_GI'                  => 'Barley TGI ESTCl.',
        'BARLEY_TUG'                 => 'Barley TUG ESTCl.',
        'BARLEY_ESTCLUSTER_PLANTGDB' => 'Barley ESTCl. PlantGDB',
        'BRASSICA_EST' => '',                   #Arabidopsis - disable
        'MAIZE_CDS'    => 'Maize CDSs',
        'MAIZE_EST'    => 'Maize ESTs',
        'MAIZE_GI'     => 'Maize TGI ESTCl.',
        'MAIZE_TUG'    => 'Maize TUG ESTCl.',
        'MAIZE_CORNSENSUS'      => 'Maize ESTCl. MMP Cornsensus',
        'MILLET_EST'            => 'Millet EST',
        'RYEGRASS_EST'          => 'Ryegrass EST Vialactia',
        'RYEGRASS_CLUSTER'      => 'Ryegrass ESTCl. Vialactia',
        'SORGHUM_EST'           => 'Sorghum EST',
        'SORGHUM_GI'            => 'Sorghum TGI ESTCl.',
        'SORGHUM_TUG'           => 'Sorghum TUG ESTCl.',
        'SORGHUM_CLUSTER_PRATT' => 'Sorghum Pratt ESTCl.',
        'SUGARCANE_EST'         => 'Sugarcane EST',
        'WHEAT_EST'             => 'Wheat EST',
        'WHEAT_GI'              => 'Wheat ESTCl. TGI',
        'WHEAT_TUG'             => 'Wheat ESTCl. TUG',

        #GSS Tracks
        'RICE_JAPONICA_BACEND'           => 'O.ja. BES IRGSP',
        'OR_BBA'                         => 'O.ni. BES OMAP',
        'OR_CBA'                         => 'O.ru. BES OMAP',
        'RICE_BRACHYANTHA_BACEND'        => 'O.br. BES OMAP',
        'RICE_ALTA_BACEND'               => 'O.al. BES OMAP',
        'RICE_AUSTRALIENSIS_BACEND'      => 'O.al. BES OMAP',
        'RICE_GLABERRIMA_BACEND'         => 'O.gl. BES OMAP',
        'RICE_NIVARA_BACEND'             => 'O.ni. BES OMAP',
        'RICE_PUNCTATA_BACEND'           => 'O.pu. BES OMAP',
        'RICE_RUFIPOGON_BACEND'          => 'O.ru. BES OMAP',
        'MAIZE_HI_COT_BENNETZEN'         => 'Maize HC Bennetzen',
        'MAIZE_METH_FILT_TIGR'           => 'Maize MF Orion',
        'MAIZE_METH_FILT_CSHL/MCCOMBIE'  => 'Maize MF CSHL',
        'MAIZE_METH_FILT_CSHL_MCCOMBIE'  => 'Maize MF CSHL',
        'MAIZE_HI_COT_TIGR'              => 'Maize MF TIGR',
        'MAIZE_METH_FILT_HI_COT_CLUSTER' => 'Maize HC-MF TIGR',
        'RYEGRASS_SEQUENCE'              => 'Ryegrass MF Orion',
        'RYEGRASS_ASSEMBLY'              => 'Ryegrass MF Cluster Orion',
        'SORGHUM_RYEGRASS_ASSEMBLY'      => 'Ryegrass MF Cluster Orion',
        'SORGHUM_RYEGRASS_SEQ'           => 'Ryegrass MF Orion',
        'SORGHUM_GSS-READ_KLEIN'         => 'Sorghum GSS Klein',
        'SORGHUM_ORION'                  => 'Sorghum MF Orion',

        # FST tracks
        'TOS17'                       => 'Rice FST Tos17',
        'RICE_TOS17_INSERT'           => 'Rice FST Tos17',
        'RICE_DS_INSERT'              => 'Rice FST Ds',
        'RICE_T_DNA_INSERT'           => 'Rice FST T-DNA',
        'RICE_TRANSPOSON_INSERT_SITE' => 'Rice FST IS',
        'MAIZE_MU_INSERT'             => 'Maize FST Mu',

        # Array tracks
        'WHEAT_CONSENSUS'             => 'Wheat Cons. A61K',
        'RICE_CONSENSUS_AFF'          => 'Rice Cons. A57K',
        'RICE_ARRAYCONSENSUS_AFFY57K' => 'Rice Cons. A57K',
        'MAIZE_CONSENSUS_AFF'         => 'Maize Cons. A18K',
        'WHEAT_CONSENSUS_AFF'         => 'Wheat Cons. A61K',
        'BARLEY1_GENECHIP_EXEMPLARS'  => 'Barley GeneChip A22K',
        'RICE_ESTOLIGO_TIGR'          => 'Rice Oligo NSF20K',
        'RICE_MPSS_BLAKE'             => 'Rice MPSS Blake',

        # Marker tracks
        'RICE_MARKER'               => 'Rice RFLP',
        'RICE_RFLP_MARKER'          => 'Rice RFLP',
        'RICE_RFLP_MARKER_RICE'     => 'Rice RFLP',
        'RICE_RFLP_MARKER_NON_RICE' => 'Non-Rice RFLP',
        'RICE_SSR'                  => 'Rice SSR',

        # MDRs
        'MDR_0.25' => 'MDR repeats=2',
        'MDR_1'    => 'MDR repeats=10',
        'MDR_2'    => 'MDR repeats=100',
        'MDR_3'    => 'MDR repeats=1000',

        # MAIZE CytoView Tracks
        'MAIZE_MARKER'    => 'Maize Markers',
        'MAIZE-OVERGOS'   => 'Maize Overgos',
        'OVERGO-AP'       => 'Overgos/Paterson',
        'OVERGO-DUPONT'   => 'Overgos/DuPont',
        'PCR'             => 'PCR Markers',
        'ELECTRONIC_SSR'  => 'Electronic SSRs',
        'CORE_BIN_MARKER' => 'Core Markers',
        'CENTA-INT7'      => 'Probes/CentA-int7',
        'P-CENTC'         => 'Probes/CentC',
        'CENTA-LTR'       => 'Probes/CentA-LTR',
        'CENT-LTR'        => 'Probes/Cent-LTR',
        'TELO'            => 'Probes/Telomeric',
        'RIBO'            => 'Probes/Ribosomal',
        'CHLORO'          => 'Probes/Chloroplasts',
        'KNOB'            => 'Probes/Knob',
        'MITO'            => 'Probes/Mitochondrial',

        'BARLEY_ARRAYCONSENSUS_AFFY22K' => 'Barley Cons. A22K',
        'BARLEY_ESTCLUSTER_PLANTGDB'    => 'Barley PlantGDB ESTCl.',
        'MAIZE_ARRAYCONSENSUS_AFFY18K'  => 'Maize Cons. A18K',
        'MAIZE_ARRAYTARGET_NSF58K'      => 'Maize Target NSF58K',
        'MAIZE_BACEND'                  => 'Maize BAC Ends',
        'MAIZE_CONSENSUS'               => 'Maize Consensus',
        'MAIZE_ESTCLUSTER_PLANTGDB'     => 'Maize PlantGDB ESTCl.',
        'MAIZE_MAGI_ISU'                => 'Maize ISU MAGI',
        'MAIZE_MARKERS'                 => 'Maize FPC Markers',
        'MAIZE_MRNA'                    => 'Maize mRNA',
        'MAIZE_SEQ_PANZEA'              => 'PanZea Sequences',
        'MAIZE_WGS_JGI'                 => 'Maize JGI WGS',
        'NONRICE_RFLP_MARKER'           => 'Non-Rice RFLP Markers',
        'OTHER-POACEAE_EST'             => 'Other Poaceae ESTs',
        'RICEALTA_BACEND_OMAP'          => 'O.al. BES OMAP',
        'RICEAUSTRALIENSIS_BACEND_OMAP' => 'O.au. BES OMAP',
        'RICE_BAC'                      => 'Rice BACs',
        'RICEBRACHYANTHA_BACEND_OMAP'   => 'O.br. BES OMAP',
        'RICECOARCTATA_BACEND_OMAP'     => 'O.co. BES OMAP',
        'RICE_ESTCLUSTER_PLANTGDB'      => 'Rice PlantGDB ESTCl.',
        'RICE_FSTTRANSPOSON'            => 'Rice FST Transposons',
        'RICEGLABERRIMA_BACEND_OMAP'    => 'O.gl. BES OMAP',
        'RICEGRANULATA_BACEND_OMAP'     => 'O.gr. BES OMAP',
        'RICEJAPONICA_BACEND_OMAP'      => 'O.ja. BES OMAP',
        'RICEMINUTA_BACEND_OMAP'        => 'O.mi. BES OMAP',
        'RICENIVARA_BACEND_OMAP'        => 'O.ni. BES OMAP',
        'RICEOFFICINALIS_BACEND_OMAP'   => 'O.of. BES OMAP',
        'RICEPUNCTATA_BACEND_OMAP'      => 'O.pu. BES OMAP',
        'RICERIDLEYI_BACEND_OMAP'       => 'O.ri. BES OMAP',
        'RICERUFIPOGON_BACEND_OMAP'     => 'O.ru. BES OMAP',
        'SORGHUM_ESTCLUSTER3P_LGBPRATT' => 'Sorghum Pratt ESTCl.',
        'SORGHUM_ESTCLUSTER_PLANTGDB'   => 'Sorghum PlantGDB ESTCl.',
        'SORGHUM_MARKERS'               => 'Sorghum Markers',
        'WHEAT_ARRAYCONSENSUS_AFFY61K'  => 'Wheat Cons. A61K',
        'WHEAT_ESTCLUSTER_PLANTGDB'     => 'Wheat PlantGDB ESTCl.',
        'WHEAT_MARKERS'                 => 'Wheat Markers',

        # Barley_est
        # Barley_GI
        # Maize_est
        # Maize_GI
        # Maize_hi_cot_Bennetzen
        # Maize_hi_cot_TIGR
        # Maize_meth_filt_CSHL_Mccombie
        # Maize_meth_filt_hi_cot_cluster
        # Maize_meth_filt_TIGR
        # Maize_Mu_Insert
        # Millet_est
        # Rice_ArrayConsensus_Affy57K
        # Rice_est
        # Rice_GI
        # Rice_ind_cluster
        # Rice_ind_est
        # Rice_jap_cDNA_KOME
        # Rice_mRNA
        # Rice_rflp_marker
        # Rice_T_DNA_Insert
        # Rice_tos17_insert
        # Ryegrass_Assembly
        # Ryegrass_Sequence
        # Sorghum_CDNA
        # Sorghum_est
        # Sorghum_GI
        # Sorghum_gss
        # Sorghum_orion
        # Sugarcane_est
        # Wheat_est
        # Wheat_GI
    );
    my $label = $labels{ uc($logic_name) };
    return defined($label) ? $label : $logic_name;
}

sub ADD_ALL_PROTVIEW_DOMAIN_TRACKS {
    my $self = shift;
    my $POS = shift || 25000;

    # Only add once!
    if ($self->{_ADDED_ALL_PROTVIEW_DOMAIN_TRACKS}) { return $POS }
    $self->{_ADDED_ALL_PROTVIEW_DOMAIN_TRACKS}++;

    #  my %logic_names;
    my @analysis_data;
    for my $sp (@{ $self->{'species_defs'}->ENSEMBL_SPECIES }) {
        my $anal_adapt = $self->registry->get_adaptor($sp, 'core', 'Analysis')
            || next;
        my $class = 'ProteinFeature';
        for my $a (@{ $anal_adapt->fetch_all_by_feature_class($class) }) {
            my $logic_name = $a->logic_name;
            push @analysis_data,
                {
                analysis   => $a,
                species    => $sp,
                logic_name => $logic_name,
                lname      => uc($logic_name),
                label      => $a->display_label || $logic_name
                };
        }
    }

    my %colours = (
        'BLASTPRODOM' => 'purple',
        'HMMPIR'      => 'red',
        'HMMSMART'    => 'darkgreen',
        'HMMTIGR'     => 'green',
        'PFAM'        => 'grey33',
        'PFSCAN'      => 'orange',
        'PRINTS'      => 'rust',
        'NCOILS'      => 'darkblue',
        'SEG'         => 'gold2',
        'SCANPROSITE' => 'orange',
        'SCANREGEXP'  => 'orange',
        'SUPERFAMILY' => 'blue',
    );

    my %no_caption = (
        'SEG'    => 1,
        'NCOILS' => 1
    );

    # Sort by logic, but send %no_caption tracks to end
    for my $a_data (
        sort {
                   ($no_caption{ $b->{lname} } <=> $no_caption{ $a->{lname} })
                || ($a->{label} cmp $b->{label})
        } @analysis_data
        )
    {

        my $analysis = $a_data->{analysis};
        my $lname    = $a_data->{lname};
        my $colour   = $colours{$lname} || 'black';

        my $track_label = ($a_data->{label});

        my $compact = $no_caption{$lname} ? 1 : 0;

        # Would like to be able to have different labels on a per-species
        # basis, but this is not possible with Ensembl's current approach!
        $self->add_new_track_Pgeneric_match($a_data->{logic_name},
            $a_data->{label}, $colour, $POS++, $compact);
    }

}

sub add_GSV_protein_domain_track {
    my ($self, $code, $text_label, $pos, %pars) = @_;
    $self->add_track(
        $code,
        'on'         => 'on',
        'pos'        => $pos,
        'glyphset'   => 'GSV_generic_domain',
        '_menu'      => 'features',
        'available'  => "features $code",
        'logic_name' => $code,
        'caption'    => $text_label,
        'dep'        => 20,
        'url_key'    => uc($code),
        'colours' => { $self->{'_colourmap'}->colourSet('protein_features') },
        %pars
    );
}

sub storable : lvalue {
    ### a
    ### Set whether this ScriptConfig is changeable by the User, and hence
    ### needs to access the database to set storable do
    ### $script_config->storable = 1; in SC code...
    $_[0]->{'storable'};
}

sub altered : lvalue {
    ### a
    ### Set to one if the configuration has been updated...
    $_[0]->{'altered'};
}

sub get_user_settings {
    my $self = shift;
    return $self->{'user'};
}

1;
