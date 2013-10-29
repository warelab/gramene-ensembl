package Maize::Page;

# $Id: Page.pm,v 1.12 2007-09-12 01:45:36 pasternak Exp $
# vim: set ts=2:

=head1 NAME

Maize::Page - Maize page wrapper

=head1 DESCRIPTION

Maize::Page creates the headers and footers and site
navigation bars around both dynamic and static content for all Maize
pages.

=head1 USAGE

Note that PerlSetEnv WRAP_ONLY implies that WrapHTML will ignore
all urls that don't include this value in the filename part of 
the path.

=head2 tutorials and frames

Tutorial = powerpoint that was converted to html by Open Office.
This html uses frames and javascript.
Each tutorial has its own directory.
Only the main .html that defines the frames is called *tutorial*
The WRAP_ONLY PerlSetEnv  is set to 'tutorial' so that WrapHTML will
alter only this file and the ___tutorial__EMPTY.html url (see below)
And not add logo & navbar to all frame components.

On pages that contain <frameset>, WrapHTML does this:
		Add an 'empty' frame at the top to contain the logo & navbar
				For /myloc/myframe.html this is /myloc/myframe.htmlEMPTY.html
					 which is actually read from /empty.html
							(its <base> makes it's link break out of frames)
		    Com
		Treat the <noframes> <body> as a normal page



=cut

use strict;
use CGI;
use Carp;

use constant CELL_LINEPADDING   => 0;
use constant CELL_LINESPACING   => 0;
use constant DEFAULT_PAGEWIDTH  => 500;
use constant DEFAULT_PANELWIDTH => 100;

# HTML Templates
use constant DIV_TEMPLATE 
    => qq(<div id="%s">
%s
</div>
);

use constant HTML_OPTION_TEMPLATE
    => qq(<option value="%s">%s</option>
);

use constant CSS_OPTION_TEMPLATE
    => qq(<link rel="stylesheet" href="%s" type="text/css"\\>
);

use constant JS_OPTION_TEMPLATE
    => qq(<script type="text/javascript" src="%s"></script>
);

use constant SEARCH_OPTION_TEMPLATE
    => qq(
<form METHOD=GET ACTION="/db/searches/browser" style="display:inline">
<select NAME="search_type">%s
</select>
<input TYPE="text" NAME="query" size="15" maxlength="100">
<input TYPE="submit" value="Search">
<input TYPE="hidden" value="Long">
<input TYPE="hidden"   name="RGN" value="on">
</form>);

use constant GUIDE_OPTION_TEMPLATE
    => qq (
<table border="0" cellpadding="0" cellspacing="0" width="100%%">
 <tr>
  <!-- Maize banner + logo -->
  <td rowspan="2">%s<br>%s%s</td>
  <!-- Ensembl logo + species -->
  <td valign="bottom" colspan="3" width="100%%">&nbsp;%s</td>
 </tr>
 <tr>
  <!-- Search -->
  <td align="right">%s</td>
  <!-- Feedback -->
  <td>%s</td>
 </tr>
</table>
);

use constant LOGO_DIV
    => qq(
<a href="/">
<div id="maize_logo">&nbsp;</div>
</a>
    );


use constant IMGLINK_OPTION_TEMPLATE
    => qq(<a href="%s"><img src="%s" alt="%s" title="%s" border="0"></a>);


my %PANELS;
my $bottom_menu_bar;

sub new {
    my $class     = shift;
    my $request   = shift || Apache2::RequestUtil->request();

    my $hostname = $request->get_server_name
        . ($request->get_server_port == 80 ? '' : ":" . $request->get_server_port);

    (my $current_url = $request->uri) =~ s,^//*,/,;
    $current_url =~ s!/perl/$ENV{ENSEMBL_SPECIES}/!/$ENV{ENSEMBL_SPECIES}/!
        if defined $ENV{ENSEMBL_SPECIES};

    $current_url =~ s!EMPTY.html$!!;

    return bless {
        r            => $request,
        footer_labels => $request->dir_config('FooterLabels'),
        hostname     => $hostname,
        current_url  => $current_url,
        current_args => scalar($request->args),
        stylesheet   => $request->dir_config('Stylesheet'),
        stylesheets  => [ split(':', $request->dir_config('Stylesheet')) ],
        javascripts  => [ split(':', $request->dir_config('JavaScript')) ],
        background   => $request->dir_config('Background'),
        bgcolor      => $request->dir_config('Bgcolor'),
        logo         => $request->dir_config('Logo'),
        enslogo      => $request->dir_config('EnsLogo'),
        banner       => $request->dir_config('Banner'),
        pagewidth    => $request->dir_config('PageWidth') || DEFAULT_PAGEWIDTH,
        panelwidth   => $request->dir_config('PanelWidth') || DEFAULT_PANELWIDTH,
        footer       => $request->dir_config('Footer'),
        sidebar      => $request->dir_config('SideBar'),
    }, $class;
} ## end sub new

sub panel      { shift->{panel} }
sub banner     { shift->{banner} }
sub logo       { shift->{logo} }
sub enslogo    { shift->{enslogo} }
sub background { shift->{background} }
sub bgcolor    { shift->{bgcolor} }

sub stylesheet {
    warn("DEPRECATED; use stylesheets instead");
    warn "From file: ", join(" line ", (caller(0))[ 1 .. 2 ]);
    shift->{stylesheet};
}
sub stylesheets { @{ shift->{stylesheets} || [] } }
sub javascripts { @{ shift->{javascripts} || [] } }
sub pagewidth   { shift->{pagewidth} }
sub panelwidth  { shift->{panelwidth} }
sub footer      { shift->{footer} }
sub sidebar     { shift->{sidebar} }
sub hostname    { shift->{hostname} }
sub current_url { shift->{current_url} }
sub current_args { shift->{current_args} }
sub modified     { shift->panel->modified }
sub r            { shift->{r} }

sub url_for_form {
    my $self         = shift;
    my $current_args = $self->current_args;
    return $self->current_url . ($current_args ? '?' . $current_args : '');
} ## end sub url_for_form

sub stylesheet_link {
    my $self      = shift;
    my @sheets    = $self->stylesheets;
    my @scripts   = $self->javascripts;
    return
        join("", map { sprintf(CSS_OPTION_TEMPLATE, $_) } @sheets)
        . join("", map { sprintf(JS_OPTION_TEMPLATE, $_) } @scripts)
        . qq(<link rel="Shortcut Icon" type="image/ico" href="/favicon.ico" />);
}

#
# Starts the HTML for the page, including stylesheets etc needed
# for the header. Is a wrapper for CGI::start_html, and you can pass
# appropriate args if wanted.
#
sub start_html {
    my $self = shift;
    my %args = @_;
    $args{-style} = [ map { { src => $_ } } $self->stylesheets() ];
    $args{-script} = [ map { { src => $_, language => 'javascript' } }
            $self->javascripts() ];
    return CGI::start_html(%args);
}

#
# This starts the table that includes the navigation panel and
# background.
#
sub start_body {
    my $self = shift;
    my %attr = @_;      # { -ensembl=> .., -bodyattr=>'..', -bodyfirst=>'...' }
    my $background = $self->background;
    my $bgcolor    = $self->bgcolor;
    my $pagewidth  = $self->pagewidth;
    my $panelwidth = $self->panelwidth;
    my $ensembl    = $attr{-ensembl} || $ENV{ENSEMBL_SPECIES} || 0;
    my $bodyattr   = $attr{-bodyattr} || '';
    my $bodyfirst  = $attr{-bodyfirst} || '';

    my $text = $bodyattr ? qq(<body $bodyattr>\n) : qq(<body>\n);
    $text .= $bodyfirst;

    my $current_url = $self->current_url;
    my $live_link   = "http://www.maizesequence.org/"
        . $current_url
        . ($self->current_args ? '?' . $self->current_args : '');

    $text .= qq(<a name="top"></a>);
    return $text;
} ## end sub start_body

sub end_body {
    my $self = shift;

    return slurp($self->footer());
}

=pod

=head2 make_sidebar
    Create the sidebar

=cut

sub make_sidebar {
    my $self = shift;
    return slurp($self->sidebar());
}

=pod

=head2 slurp
    Returns the contents of the file

=cut

sub slurp {
    my ($filename) = @_;
    my $output;
    local $/;
    open my $file, '<', $filename or croak "Cannot open $filename: $!\n";  
    $output = <$file>;
    close $file or croak "Cannot close filehandle for $filename: $!\n";
    return $output;
}

=pod

=head2 create_search_form
    Creates a Search HTML form

=cut

sub create_search_form {
    my $self = shift;
    # ---Search---

    my @options_html = (sprintf(HTML_OPTION_TEMPLATE, 'All', 'Find anything'));
    my @options = qw(Markers);
    for my $option (@options) {
        push @options_html, sprintf(HTML_OPTION_TEMPLATE, $option, $option);
    }

    return sprintf(SEARCH_OPTION_TEMPLATE, join('', @options_html));
}

=pod

=head2 title
    Returns the title of this page

=cut

sub title {
    my $self = shift;
    return "TITLE";
}

1;

=pod

=head1 AUTHOR

Lincoln Stein.  Modified by a cast of thousands.

=cut

