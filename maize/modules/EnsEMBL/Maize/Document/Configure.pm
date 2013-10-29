package EnsEMBL::Maize::Document::Configure;

use Maize::Page;
use EnsEMBL::Web::Form;
use EnsEMBL::Maize::Document::HTML::GenomeEntryPoints;

use EnsEMBL::Web::Root;
our @ISA = qw(EnsEMBL::Web::Root);

sub common_page_elements {
    my ($self, $document) = @_;

    # We have to attach some variables to the apache request as these are
    # required by Maize::Page. The variables themselves are set in the
    # SiteDefs.pm plugin.
    # These should remain available on the Apache Request within other modules
    my $r = Apache2::RequestUtil->request;
    foreach my $key (keys %SiteDefs::PerlSetVar) {
        $r->dir_config->add($key, $SiteDefs::PerlSetVar{$key});
    }

    # Maize header
    $document->replace_body_element(
        'masthead' => 'EnsEMBL::Maize::Document::HTML::MastHead');

    $document->replace_body_element(
        'helplink' => 'EnsEMBL::Maize::Document::HTML::HelpLink');

    $document->replace_body_element(
        'searchbox' => 'EnsEMBL::Maize::Document::HTML::SearchBox');

    # Maize footer
    $document->replace_body_element(
        'copyright' => 'EnsEMBL::Maize::Document::HTML::Copyright');

    $document->replace_body_element(
        'release' => 'EnsEMBL::Maize::Document::HTML::Null');

    return 1;
}

sub static_page_elements {
    my $self = shift;
    my ($document) = @_;

    # Maize content
    $document->replace_body_element(
        'content' => 'EnsEMBL::Maize::Document::HTML::Content');

    return 1;
}

sub extra_configuration {
    my ($self, $document) = @_;

    # MAIZE SPECIFIC stylesheets and javascript files
    my $page = Maize::Page->new()
        || die("Cannot create a Maize::Page");
    foreach my $css ($page->stylesheets) {
        $document->stylesheet->add_sheet('screen', $css);
    }
    foreach my $file ($page->javascripts) {
        $document->javascript->add_source($file);
    }
}

sub common_menu_items {
    my ($self, $document) = @_;
    my $menu = $document->menu;

    # Wipe all existing blocks
    map { $menu->delete_block($_) } $menu->blocks;

    $menu->add_block('info', 'attribute', 'Site Info');
    $menu->add_entry(
        'info',
        'code'  => 'home',
        'href'  => '/',
        'text'  => 'Home',
        'title' => 'Go to browser home',
    );
    $menu->add_entry(
        'info',
        'code'  => 'overview',
        'href'  => '/overview.html',
        'text'  => 'Overview',
        'title' => 'Project Overview',
    );
    $menu->add_entry(
        'info',
        'code'  => 'faq',
        'href'  => '/faq.html',
        'text'  => 'FAQ',
        'title' => 'Frequenty Asked Questions',
    );
    $menu->add_entry(
        'info',
        'code'  => 'feedback',
        'href'  => '/perl/feedback',
        'text'  => 'Feedback',
        'title' => 'Provide Feedback',
    );
    $menu->add_entry(
        'info',
        'code'  => 'blastview',
        'href'  => '/Multi/blastview',
        'text'  => 'BLAST',
        'title' => 'Run a Sequence Similarity (BLAST) Search',
    );
    $menu->add_entry(
        'info',
        'code'  => 'ftp',
        'href'  => 'http://ftp.maizesequence.org',
        'text'  => 'FTP',
        'title' => 'FTP downloads of all maize-related data',
        'raw'   => 1,
    );
    $menu->add_entry(
        'info',
        'code'  => 'notification',
        'href'  => '/perl/notification',
        'text'  => 'BAC Notification',
        'title' => 'Subscribe to data updates via RSS feeds',
        'raw'   => 1,
    );

    my $gramene_url = $SiteDefs::PerlSetVar{GrameneUrl};
    my $options     = [];
    for my $genome (@{ $SiteDefs::PerlSetVar{GrameneGenomes} }) {
        push @$options,
            +{
            'text' => $genome->[0],
            'href' => "$gramene_url/$genome->[1]",
            };
    }
    $menu->add_entry(
        'info',
        'code'    => 'gramene',
        'href'    => $gramene_url,
        'text'    => 'Gramene',
        'title'   => 'Jump to another genome on Gramene',
        'options' => $options,
    );

    $self->add_genome_navigation($document);
}

sub dynamic_menu_items {
    my $self = shift;
    my ($document) = @_;
}

=pod

=head2 add_genome_navigation
    Adds click-thrus for clones and contigs from database

=cut

sub add_genome_navigation {
    my $self = shift;
    my ($document) = @_;

    use CGI;
    my $cgi = new CGI;
    if ($cgi->param('dev')) {
        return;
    }
    my $panel = EnsEMBL::Maize::Document::HTML::GenomeEntryPoints->new(
        { 'all' => 1 });

    my $menu = $document->menu;
    $menu->add_block('navigation', 'raw', 'Navigation',
        'html' => $panel->render);
    $menu->change_block_attribute('navigation', 'priority', 5);
}

1;
