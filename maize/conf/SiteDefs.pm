package EnsEMBL::Maize::SiteDefs;
use strict;

# These are the maize-specific edits to the main Ensembl SiteDefs.pm file
sub update_conf {
    warn("===> MAIZE");
    $SiteDefs::ENSEMBL_PORT = 8081;

# $SiteDefs::ENSEMBL_PROXY_PORT = 80;   # Port used for self-referential URLs:

    $SiteDefs::ENSEMBL_USER  = getpwuid($>);
    $SiteDefs::ENSEMBL_GROUP = getgrgid($));

    $SiteDefs::ENSEMBL_SERVERADMIN = 'webmaster@maizesequence.org';
    $SiteDefs::ENSEMBL_SERVERNAME  = 'dev.maizesequence.org';
    $SiteDefs::ENSEMBL_MAIL_SERVER = 'maizesequence.org';
    $SiteDefs::ENSEMBL_MAIL_ERRORS = 1;
    $SiteDefs::ENSEMBL_ERRORS_TO   = 'webmaster@maizesequence.org';

    $SiteDefs::ENSEMBL_DEBUG_FLAGS = 24;

    $SiteDefs::ENSEMBL_LONGPROCESS_MINTIME = 10;

    $SiteDefs::ENSEMBL_USERDB_NAME = 'ensembl_web_user_db_46';
    
    $SiteDefs::ENSEMBL_LOGINS = 1;

    #----------
    # Logging
    my $LOG_ROOT = $SiteDefs::ENSEMBL_SERVERROOT . "/logs";
    $SiteDefs::ENSEMBL_PIDFILE   = "$LOG_ROOT/httpd.pid";
    $SiteDefs::ENSEMBL_ERRORLOG  = "$LOG_ROOT/error.log";
    $SiteDefs::ENSEMBL_CUSTOMLOG = "$LOG_ROOT/access.log combined";

    #----------
    # Species stuff
    $SiteDefs::ENSEMBL_PERL_SPECIES      = 'Zea_mays2';    # Default species
    $SiteDefs::ENSEMBL_PRIMARY_SPECIES   = 'Zea_mays2';     # Default species
    $SiteDefs::ENSEMBL_SECONDARY_SPECIES = 'Zea_mays';    # Default species
    %SiteDefs::__species_aliases         = (
        'Zea_mays'          => [ ('maize',          'zm') ],
        'Zea_mays2'         => [ ('Zea_mays_clone', 'zmays', 'zm2') ],
        'Zea_mays_external' => [ ('maize_old',      'zm3') ],
        'Oryza_sativa_japonica' =>
            [ ('rice', 'osativa', 'os', 'Oryza_sativa') ],
    );

    $SiteDefs::ENSEMBL_MART_ENABLED = 0;    #Turn off mart for now

    my $MAIZE_ROOT     = get_maize_root($SiteDefs::ENSEMBL_PLUGINS);
    my $MAIZE_DOCSROOT = "$MAIZE_ROOT/htdocs";

    unshift(@SiteDefs::ENSEMBL_LIB_DIRS, "$MAIZE_ROOT/modules");

    @SiteDefs::ENSEMBL_LIB_DIRS
        = grep { !m#gramene-live/lib/perl# } @SiteDefs::ENSEMBL_LIB_DIRS;
    $SiteDefs::PerlSetVar{Stylesheet} = join(':',
        map {"/stylesheets/$_.css"} ('maize', 'maize-page', 'nifty'));
    $SiteDefs::PerlSetVar{JavaScript} = join(':',
        map {"/js/$_.js"}
            ('nifty', 'forms', 'drag_imagemap', 'prototype', 'effects'));
    $SiteDefs::PerlSetVar{Footer}     = "$MAIZE_DOCSROOT/layout/footer.html";
    $SiteDefs::PerlSetVar{SideBar}    = "$MAIZE_DOCSROOT/layout/sidebar.html";
    $SiteDefs::PerlSetVar{Background} = q{};
    $SiteDefs::PerlSetVar{Bgcolor}    = q{};
    $SiteDefs::PerlSetVar{Logo}       = q{};
    $SiteDefs::PerlSetVar{EnsLogo}    = q{};
    $SiteDefs::PerlSetVar{Banner}     = q{};
    $SiteDefs::PerlSetVar{PageWidth}  = q{};
    $SiteDefs::PerlSetVar{PanelWidth} = q{};

    $SiteDefs::PerlSetVar{GrameneUrl}     = 'http://dev.gramene.org';
    $SiteDefs::PerlSetVar{GrameneGenomes} = [
        [ 'Oryza sativa'         => 'Oryza_sativa' ],
        [ 'Oryza rufipogon'      => 'Oryza_rufipogon' ],
        [ 'Arabidopsis thaliana' => 'Arabidopsis_thaliana' ],
        [ 'Organelles'           => 'organelles' ],
    ];
}

=pod

=head2 get_maize_root
    Extract maize root

=cut

sub get_maize_root {
    my ($plugins) = @_;
    my @temporary_array = @$plugins;
    while (my ($plugin, $root) = splice(@temporary_array, 0, 2)) {
        if ($plugin eq 'EnsEMBL::Maize') {
            return $root;
        }
    }
    die("Unable to locate Maize root!\n");
}

1;
