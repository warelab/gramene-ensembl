package EnsEMBL::HTML::Page;

#
# EnsEMBL HTML page headers and footers
#

use strict;
use Apache;
use Apache::Constants qw(:response :methods :http);
use Apache 'exit';
use Sys::Hostname;
use Exporter;

#use Sanger::Graphics::JSTools;

# Maize page header
use Maize::Page;

use vars qw( @ISA @EXPORT $SPECIES_DEFS);

use Exporter();
@ISA = qw(Exporter);

# Default Exports
@EXPORT = (
    'ensembl_exit',
    'ensembl_exception',
    'ensembl_warning',
    'ensembl_http_header',     # Renamed from http_header
    'ensembl_page_header',     # Renamed from make_cgi_header
    'ensembl_page_footer',     # Renamed from make_cgi_footer
    'ensembl_search_table',    # Renamed from print_form
    'ensembl_help_link',       # Renamed from helpLink
    'ensembl_affiliation_table',
    'ensembl_navigation_table',
);

#BEGIN:{
#  $SPECIES_DEFS = SpeciesDefs->new();
#  $SPECIES_DEFS->retrieve() ||
#    &ensembl_exception("Unable to retrieve species defs");
#}

#----------------------------------------------------------------------
# Creates a new EnsEMBL::HTML::Page object
sub new {
    my $class = shift;
    my $self  = bless {}, $class;

    $self->species($ENV{ENSEMBL_SPECIES});
    $self->script($ENV{ENSEMBL_SCRIPT});
    $self->{_r} = Apache->request;

    return $self;
}

#----------------------------------------------------------------------
# Flush apache output buffer. Wrapper for Apache request rflush method
sub flush {
    my $self = shift;
    $self->{_r}->rflush;
}

#----------------------------------------------------------------------
# Species accessor
sub species {
    my $key  = '_species';
    my $self = shift;
    if (@_) {
        my $species = shift;
        $self->{$key} = $species;
        $self->species_common_name($species);
    }
    return $self->{$key};
}

#----------------------------------------------------------------------
# Species common name accessor
sub species_common_name {
    return '';

    #  my $key = '_species_common_name';
    #  my $self = shift;
    #  if( @_ ){
    #    my $species = shift;
    #    local $SIG{__WARN__} = sub{}; # Trap conf-not-found errors
    #    my $species_common = $SPECIES_DEFS->get_config( $species,
    #						    'SPECIES_COMMON_NAME' );
    #    $self->{$key} = $species_common || $species;
    #  }
    #  return $self->{$key};
}

#----------------------------------------------------------------------
# Script accessor
sub script {
    my $key  = '_script';
    my $self = shift;
    if (@_) {
        my $script = shift;
        $self->{$key} = lc($script);
        $self->script_pretty_name($script);
    }
    return $self->{$key};
}

#----------------------------------------------------------------------
# Script pretty name accessor
sub script_pretty_name {
    my $key  = '_script_pretty_name';
    my $self = shift;
    if (@_) {
        my $name = lc(shift @_);
        $name = ucfirst $name;
        $name =~ s/view/View/;
        $self->{$key} = $name;
    }
    return $self->{$key};
}

#----------------------------------------------------------------------
# Adds javascript file to the page header
sub add_javascript_file {
    my $self     = shift;
    my $filename = shift;
    $self->{_javascript_files} ||= [];
    push @{ $self->{_javascript_files} }, $filename;
}

#----------------------------------------------------------------------
# Adds jscript code to the page header
sub add_javascript_to_header {
    my $self = shift;
    my $code = shift || return;
    $self->{_javascript_in_header} ||= [];
    push @{ $self->{_javascript_in_header} }, $code;
}

#----------------------------------------------------------------------
# Adds javascript to the body.
sub add_javascript_to_body {
    my $self = shift;
    my $code = shift || return;
    $self->{_javascript_in_body} ||= [];
    push @{ $self->{_javascript_in_body} }, $code;
}

#----------------------------------------------------------------------
# Adds a javascript onload function to the page <BODY> tag
sub add_onload_function {
    my $self = shift;
    my $function = shift || return;
    $self->{_onload_functions} ||= [];
    push @{ $self->{_onload_functions} }, $function;
}

#----------------------------------------------------------------------
# Configures page to load zmenus
sub add_zmenu_javascript {
    my $self = shift;

 #  $self->add_javascript_to_header(&Sanger::Graphics::JSTools::js_menu_header);
 #  $self->add_javascript_to_body(&Sanger::Graphics::JSTools::js_menu_div);
}

#----------------------------------------------------------------------
# rapper for CGI::Header
sub http_header {
    return &CGI::header();
}

#----------------------------------------------------------------------
# Generates the HTML HEAD and BODY tags to start the page
sub start_html {
    my $self = shift;

    #--- Start HTML templates ---

    my $HTML_HEAD = qq(
<!DOCTYPE HTML PUBLIC "-//IETF//DTD HTML//EN">
<HTML>
<HEAD>
<script language="javascript">
 var LOADED = 0;
</script>
 %s
 %s
 %s
 %s
 %s
</HEAD>
<BODY %s text="#000000" bgcolor="#FFFFFF" topmargin="2" marginheight="2">
%s);

    my $HTML_TITLE = qq(
<TITLE>%s %s (%s)</TITLE>);

    my $LINK_REL = qq(<LINK rel="stylesheet" href="/gfx/EnsMart/%s">);

    my $JAVASCRIPT_FILE = qq(
  <script language="Javascript" src="%s"></script> );

    my $ONLOAD_FUNCTION = qq(onLoad="%s");

    #--- End HTML templates

    #--- Start generating HTML ---

    # Generate Page title
    my $authority = 'Maize';    # $SPECIES_DEFS->AUTHORITY || 'Ensembl';

    my $title = '';        #sprintf( $HTML_TITLE,
                           #      $SPECIES_DEFS->AUTHORITY || 'Ensembl',
                           #      $self->species_common_name,
                           #      $self->script_pretty_name, );

    # Generate stylesheet
    my $css = "EnsEMBL.css";
    $css = "" if (defined $self->{_nocss});
    $css = sprintf($LINK_REL, $css);

    # Maize page
    my $grmPAGE = Maize::Page->new(Apache->request);
    $css .= "\n" . $grmPAGE->stylesheet_link;

    # Generate base href
    my $base = "";
    if (defined $self->{_basehref}) {
        $base = qq(<base href="$self->{_basehref}" />);
    }

    #Build javascript files
    my $javascript_files = '';
    if (ref($self->{_javascript_files}) eq 'ARRAY') {
        $javascript_files = join("",
            map { sprintf($JAVASCRIPT_FILE, $_) }
                @{ $self->{_javascript_files} });
    }

    #Build header javascript
    my $header_javascript;
    if (ref($self->{_javascript_in_header}) eq 'ARRAY') {
        $header_javascript = join("\n", @{ $self->{_javascript_in_header} });
    }

    #Build body javascript
    my $body_javascript;
    if (ref($self->{_javascript_in_body}) eq 'ARRAY') {
        $body_javascript = join("\n", @{ $self->{_javascript_in_body} });
    }

    #Build onload tag attribute
    my $onload_functions = '';
    $onload_functions = sprintf($ONLOAD_FUNCTION,
        join(";", 'LOADED=1', @{ $self->{_onload_functions} || [] }));

    # Build and return the header
    return sprintf($HTML_HEAD,
        $title, $css, $base, $javascript_files, $header_javascript,
        $onload_functions, $body_javascript);
}

#----------------------------------------------------------------------
# Wrapper for exit
sub ensembl_exit {

    #use Bio::EnsEMBL::Registry; Bio::EnsEMBL::Registry->disconnect_all();
    exit;
}

#----------------------------------------------------------------------
# Prints Ensembl Exception page
#   Takes:
#	error message to display to user,
#	"real" error - i.e $@ or $! or whatever,
#	nomail - if this has a value, no email is sent.  Use for trivial "you
#	have forgotten to supply an argument" type errors.
#
sub ensembl_exception {
    my ($user_error, $exception, $nomail) = @_;

    if ($ENV{'GATEWAY_INTERFACE'} =~ /^CGI-Perl/) {

        # mod_perl: delegate error handling
        my $hostname = Sys::Hostname::hostname();
        my $r        = Apache->request();
        $r->subprocess_env->{'error_message'} = $user_error;
        $r->subprocess_env->{'error_notes'}   = $exception;

        #  $r->err_header_out('ensembl_headers_out' => 1           );
        #  $r->err_header_out('Ensembl-Error'       => $user_error );
        #  $r->err_header_out('Ensembl-Exception'   => $exception  );
        #  $r->err_header_out('Ensembl-Server'      => $hostname   );
        #  $r->err_header_out('Ensembl-Nomail'      => $nomail     );
        $r->internal_redirect('/Crash');    # Redirect to ensembl error handler
    } else {
        warn "plain exception ..>";

        # This is vanilla CGI - just barf
        print qq( 
    <h1>Ensembl Server Error</h1>
    <p>Sorry, an error occured while the server was processing your 
       request. Please email a report , quoting any additional information 
       given below, along with the URL, to Ensembl server administrator 
       using the link below</p>
    <p><b>The error was:</b></p>
    <blockquote class="error"><b><pre>$exception</pre></b></blockquote>
    <p><b>URL:</b></p>
    <blokcquote class="error"><b>$ENV{'REDIRECT_URL'}</b></blockquote>
    <p><b>HTTP Status Code:</b></p>
    <blockquote class="error"><b>$ENV{'REDIRECT_STATUS'}</b></blockquote>
    <p><b>Request Method:</b></p>
    <blockquote class="error"><b>$ENV{'REDIRECT_REQUEST_METHOD'}</b></blockquote> 
    <p><b>Query String (if known):</b></p>
    <blockquote class="error"><b>$ENV{'REDIRECT_QUERY_STRING'}</b></blockquote>
    <p><b>Error (if known):</b></p>
    <blockquote class="error"><b>$ENV{'REDIRECT_ERROR_NOTES'}</b></blockquote>
	      );
        print &ensembl_page_footer;
    }
    &ensembl_exit;
}

#----------------------------------------------------------------------
# Prints an Ensembl Warning page
sub ensembl_warning {
    my ($title, $message) = @_;
    print qq( <h3>$title</h3>\n<p>$message</p> );
    print &ensembl_page_footer;
    &ensembl_exit;
}

#----------------------------------------------------------------------
# Ensembl wrapper for CGI::Header
sub ensembl_http_header {
    return &CGI::header();
}

#----------------------------------------------------------------------
# Generates and returns headers for ensembl page
my $VIEW_IMAGE_TEMPLATE   = '/gfx/%s.gif';
my $DEFAULT_SPECIES_IMAGE = '/gfx/header/ensembl-header1.gif';

sub ensembl_page_header {

    # Get the args
    my (%opts_ref) = @_;

    my $ensembl_species = $opts_ref{ensembl_species}
        || '';    #$ENV{ENSMBL_SPECIES};
    my $common_name = $opts_ref{species_common_name} ||

        #$SPECIES_DEFS->SPECIES_COMMON_NAME ||
        $ensembl_species;

    #--- Start HTML templates ---

    my $HTML_HEAD = qq(
<!DOCTYPE HTML PUBLIC "-//IETF//DTD HTML//EN">
<HTML>
<HEAD><script language="javascript"> var LOADED = 0; </script>
 %s
 %s
 %s
 %s
</HEAD>
<BODY %s text="#000000" bgcolor="#FFFFFF" topmargin="2" marginheight="2">
%s);

    my $HTML_TITLE = qq(
<TITLE>%s %s Genome Browser (%s)</TITLE>);

    my $LINK_REL = qq(<LINK rel="stylesheet" href="/gfx/EnsMart/%s">);

    #--- End HTML templates

    #--- Start Content definition ---

    my $script_name = $ENV{'ENSEMBL_SCRIPT'};
    $script_name =~ s/.*\/(.*?)$/$1/;
    $script_name = ucfirst($script_name);
    $script_name =~ s/view$/View/;

    #--- Start javascript ---

    # Compile script
    my $script  = "";    # Java script
    my $onload  = '';
    my @onload  = ();
    my $js_divs = "";

    if (defined $opts_ref{onload}) {
        push @onload,
            (
            ref($opts_ref{'onload'}) eq 'ARRAY'
            ? @{ $opts_ref{'onload'} }
            : ($opts_ref{'onload'})
            );
    }
    push @onload, 'LOADED = 1';

#  if(defined $opts_ref{'initfocus'}) {
#    $script .=  &Sanger::Graphics::JSTools::js_init;
#    push @onload, 'init()';
#  } if(defined $opts_ref{'menus'}) {
#    $script  .= &Sanger::Graphics::JSTools::js_menu_header();
#    $js_divs .= &Sanger::Graphics::JSTools::js_menu_div();
#  } if(defined $opts_ref{'tooltips'}) {
#    $script  .= &Sanger::Graphics::JSTools::js_tooltip_header();
#    $js_divs .= &Sanger::Graphics::JSTools::js_tooltip_div();
#  } if(defined $opts_ref{'dd_menus'}) {
#    $script .= qq(\n<script language="Javascript" src="/js/dd_menus.js"></script>);
#  }
    if (defined $opts_ref{'additional_functions'}) {
        my @funcs = (
            ref($opts_ref{'additional_functions'}) eq 'ARRAY'
            ? @{ $opts_ref{'additional_functions'} }
            : ($opts_ref{'additional_functions'})
        );
        $script .= join "\n", @funcs;
    }
    if (@onload) {
        $onload = 'onLoad="' . join(';', @onload) . '"';
    }

    #--- End javascript ---

    #--- Start generating HTML ---

    my $authority = 'Maize';    # $SPECIES_DEFS->AUTHORITY;

    # Generate Page title
    my $html_title
        = sprintf($HTML_TITLE, $authority, $common_name, $script_name);

    # Generate CSS link
    my $css = "EnsEMBL.css";
    $css = "" if (defined $opts_ref{'nocss'});
    my $link_rel = sprintf($LINK_REL, $css);

    # Maize page
    my $grmPAGE = Maize::Page->new(Apache->request);
    $link_rel .= "\n" . $grmPAGE->stylesheet_link;

    # Generate base href
    my $base = "";
    if (defined $opts_ref{'basehref'}) {
        $base = qq(<base href="$opts_ref{'basehref'}" />);
    }

    # Generate the HTML header
    my $html_head = sprintf($HTML_HEAD,
        $html_title, $link_rel, $base, $script, $onload, $js_divs);

    # Generate the affiliations table
    my $affiliation_html = ensembl_affiliation_table(%opts_ref);

    # Generate the navigation table
    my $navigation_html = ensembl_navigation_table(%opts_ref);

    # Maize page
    my $grmPAGE = Maize::Page->new(Apache->request);

    # Return the bits
    return $html_head . $grmPAGE->stylesheet_link . $grmPAGE->start_body;

}

#----------------------------------------------------------------------
# Generates the page header affiliation table
sub ensembl_affiliation_table {
    return '';    # Maize - no affiliation
                  #---
                  # Get the args
    my (%opts_ref) = @_;

    my $ensembl_species = $opts_ref{ensembl_species} || $ENV{ENSEMBL_SPECIES};
    my $common_name = $opts_ref{species_common_name} ||

        #$SPECIES_DEFS->SPECIES_COMMON_NAME ||
        $ensembl_species;

    #---
    # Define templates

    use constant AFFILIATION_HTML => qq(
<TABLE width="100%%" border="0" cellspacing="0" cellpadding="0">
 <TR valign="bottom">
  <TD align="left">
    %s %s %s
  </TD>
  <TD align="right">
    %s %s
  </TD>				      
 </TR>
</TABLE> );

    my $IMAGE_HTML
        = qq(<a href="%s"><img src="%s" border="0" height="%s" alt="%s"></a>);

    my $affiliation_table_height = 40;

    #---
    # e! image

    my %ensembl_image = ();

    #( src    => $SPECIES_DEFS->AUTHORITY_LOGO,
    # height => $SPECIES_DEFS->AUTHORITY_LOGO_HEIGHT,
    # width  => $SPECIES_DEFS->AUTHORITY_LOGO_WIDTH,
    # alt    => $SPECIES_DEFS->AUTHORITY_LOGO_ALT,
    # href   => $SPECIES_DEFS->AUTHORITY_LOGO_HREF);

    #---
    # species image

    my %species_image = ();

    #$SPECIES_DEFS->HEADER_LINKS ?
    #  ( src    => $SPECIES_DEFS->HEADER_LINKS->{IMAGE1_SRC}   ||
    #              $DEFAULT_SPECIES_IMAGE,
    #    height => $SPECIES_DEFS->HEADER_LINKS->{IMAGE1_HEIGHT}||
    #              $affiliation_table_height,
    #    width  => $SPECIES_DEFS->HEADER_LINKS->{IMAGE1_WIDTH} || 174,
    #    alt    => "Ensembl $common_name",
    #    href   => $SPECIES_DEFS->HEADER_LINKS->{IMAGE1_URL}   ||
    #              "/$ensembl_species" ) : ();

    if ($opts_ref{species_image} and ref $opts_ref{species_image} eq 'HASH') {

        # User-specified values
        %species_image = %{ $opts_ref{species_image} };
    }

    #---
    # <Script>View image

    my $script = $ENV{'ENSEMBL_SCRIPT'};
    $script =~ s/.*\/(.*?)$/$1/;
    my $script_name = ucfirst($script);
    $script_name =~ s/view$/View/;
    $script = 'blank' if $script !~ /.*view$/;

    my %view_image = (
        src    => sprintf($VIEW_IMAGE_TEMPLATE, $script),
        height => $affiliation_table_height,
        width  => "",
        alt    => $script,
        href   => ''
    );

    #---
    # Sanger image

    my %sanger_image = ();

    #( src    => $SPECIES_DEFS->SANGER_LOGO ||
    #            '/gfx/header/wtsi-header.gif',
    #  height => $SPECIES_DEFS->SANGER_LOGO_HEIGHT ||
    #            $affiliation_table_height,
    #  width  => $SPECIES_DEFS->SANGER_LOGO_WIDTH ||
    #            160,
    #  alt    => $SPECIES_DEFS->SANGER_LOGO_ALT ||
    #            'The Wellcome Trust Sanger Institute',
    #  href   => $SPECIES_DEFS->SANGER_LOGO_HREF ||
    #            'http://www.sanger.ac.uk/' );

    #---
    # Collaborate image

    my %collaborate_image = ();

    #( src    => $SPECIES_DEFS->COLLABORATE_LOGO ||
    #            '/gfx/header/ebi-header.gif',
    #  height => $SPECIES_DEFS->COLLABORATE_LOGO_HEIGHT ||
    #            $affiliation_table_height,
    #  width  => $SPECIES_DEFS->COLLABORATE_LOGO_WIDTH  ||
    #            77,
    #  alt    => $SPECIES_DEFS->COLLABORATE_LOGO_ALT    ||
    #            'The European Bioinformatics Institute',
    #  href   => $SPECIES_DEFS->COLLABORATE_LOGO_HREF   ||
    #            'http://www.ebi.ac.uk/' );

    #---
    # Build and return the HTML
    my @images_html;
    foreach (
        \%ensembl_image, \%species_image, \%view_image,
        \%sanger_image,  \%collaborate_image
        )
    {
        push(
            @images_html,
            sprintf($IMAGE_HTML,
                $_->{href}, $_->{src}, $_->{height}, $_->{alt},)
        );
    }

    return sprintf(AFFILIATION_HTML, @images_html);
}

#----------------------------------------------------------------------
# Generates and returns the ensembl page header navigation table
sub ensembl_navigation_table {

    # Return Maize page page header
    my $grmPAGE = Maize::Page->new(Apache->request);
    return $grmPAGE->start_body;

    # The following is never executed for Maize
    # Get the args
    my (%opts_ref) = @_;

    my $ensembl_species = $opts_ref{ensembl_species} || $ENV{ENSEMBL_SPECIES};
    my $ensembl_link    = $opts_ref{ensembl_link}    || "/$ensembl_species/";
    my $common_name = $opts_ref{species_common_name}
        || $SPECIES_DEFS->SPECIES_COMMON_NAME
        || $ensembl_species;

    my $NAV_TABLE = qq(
<table border="0" cellspacing="0" cellpadding="0" width="100%%">
<tr><td class="header1">
<table border="0" cellspacing="0" cellpadding="0">
 <tr>
  %s
  %s
 </tr>
</table>
</td></tr></table>
);

    my $NAV_LINK = qq(<a class="trailbar" href="%s"><small>%s</small></a>);

    my $NAV_LEFT = qq(
  <td class="header1" align="center">&nbsp; %s &nbsp;</td>
  <td class="header1" align="center"><img src="/gfx/bullet.red.gif"></td>
  <td class="header1" align="center">&nbsp; %s &nbsp;</td>
  <td nowrap="nowrap" class="header2"><img src="/gfx/arrow.small.red.up.gif"></td>
 );

    my $VEGA_NAV_LEFT = qq(
   <td class="header1" align="center">&nbsp; %s &nbsp;</td> );

    my $VEGA_SPECIES_LEFT = qq(
   <td class="header1" align="center"><img src="/gfx/EnsMart/blank.gif"></td>
   <td class="header1" align="center">&nbsp; %s &nbsp;</td> );

    my $VEGA_FOCUS_SPECIES_LEFT = qq(
				   <td class="header1" align="center"><img src="/gfx/bullet.red.gif"></td>
				   <td class="header1" align="center"><b>&nbsp; %s &nbsp;</b></td>
);

    my $NAV_ENTRY = qq(
  <td nowrap="nowrap" class="header2">%s </td>
  <td nowrap="nowrap" class="header2"><img src="/gfx/arrow.small.red.up.gif"></td>
);

    # Ensembl links
    my $ens_links = {
        HOME => {
            TEXT => 'Home',
            URL  => '/'
        },
        SPECIES => {
            TEXT => "$common_name",
            URL  => "$ensembl_link"
        }
    };

    # Retrieve and parse the HEADER_LINKS config info
    my $HEADER_LINKS = undef;    #$SPECIES_DEFS->HEADER_LINKS;
    if (ref($HEADER_LINKS) ne 'HASH') { return '' }

    my %ens_order;
    foreach my $param (keys(%$HEADER_LINKS)) {
        if (my ($name, $type) = ($param =~ /^(LINK\d+)_(TEXT|URL)/o)) {
            my ($num) = ($name =~ /(\d+)/o);
            $ens_links->{$name}->{$type} = $HEADER_LINKS->{$param};
            $ens_order{$num} = $name;
        }
    }

    # Which ens links to have in the ensembl entries array
    my @ens_entries
        = map { $ens_order{$_} } sort { $a <=> $b } keys(%ens_order);
    if (!scalar(@ens_entries)) { return '' }

    # generate Vega help link
    my $ENSEMBL_SCRIPT = lc($ENV{'ENSEMBL_SCRIPT'});
    my $help_link
        = qq£javascript:void(window.open('/${ensembl_species}/helpview?se=1&kw=$ENSEMBL_SCRIPT',£;
    $help_link .= qq£'helpview','width=400,height=500,resizable,scrollbars'));£;

    # Generate the ensembl links table
    my $ensembl_html_left;

#  if ($SPECIES_DEFS->SITE_TYPE eq 'Vega') {
#      $ensembl_html_left = sprintf( $VEGA_NAV_LEFT,
#				    sprintf( $NAV_LINK,
#					     $ens_links->{HOME}->{URL},
#					     $ens_links->{HOME}->{TEXT} )
#				  );
#
#      foreach my $k (sort { $all_species{$a} cmp $all_species{$b} } keys %all_species) {
#	  if ($k eq $ensembl_species) {
#	      $ensembl_html_left .= sprintf( $VEGA_FOCUS_SPECIES_LEFT,
#					     sprintf( $NAV_LINK,
#						      "/$k/",
#						      $all_species{$k})
#					   );
#	  } else {
#	      $ensembl_html_left .= sprintf( $VEGA_SPECIES_LEFT,
#					     sprintf( $NAV_LINK,
#						      "/$k/",
#						      $all_species{$k})
#					   );
#	  }
#      }
#
#      $ensembl_html_left .= qq(
#			       <td nowrap="nowrap" class="header2"><img src="/gfx/arrow.small.red.up.gif"></td>);
#
#  } else {
    $ensembl_html_left = sprintf(
        $NAV_LEFT,
        sprintf($NAV_LINK,
            $ens_links->{HOME}->{URL},
            $ens_links->{HOME}->{TEXT}),
        sprintf($NAV_LINK,
            $ens_links->{SPECIES}->{URL},
            $ens_links->{SPECIES}->{TEXT})
    );

    #  }

    my $ensembl_html_right = '';
    foreach (@ens_entries) {
        $ensembl_html_right .= sprintf(
            $NAV_ENTRY,
            sprintf($NAV_LINK,
                $ens_links->{$_}->{URL},
                $ens_links->{$_}->{TEXT})
        );
    }

    #add vega help link to nav table
    #  if ($SPECIES_DEFS->SITE_TYPE eq 'Vega') {
    #      $ensembl_html_right .= sprintf( $NAV_ENTRY,
    #				      sprintf( $NAV_LINK,
    #					       $help_link,
    #					       "Help" ) );
    #  }
    return sprintf($NAV_TABLE, $ensembl_html_left, $ensembl_html_right);
}

#----------------------------------------------------------------------
# Generates and returns the ensembl page footer
sub ensembl_page_footer {

    # Maize-specific
    return Maize::Page->new(Apache->request)->end_body;
}

#----------------------------------------------------------------------
# Pretty-print date/time
sub format_date {
    my ($sec, $min, $hour, $day, $month, $year) = localtime();
    $year  += 1900;
    $month += 1;
    return sprintf("%04d-%02d-%02d %02d:%02d:%02d",
        $year, $month, $day, $hour, $min, $sec);
}

#----------------------------------------------------------------------
# Generates and returns the ensembl (or Vega) search table
sub ensembl_search_table {

    # Get the args
    my $q    = shift @_ || '';
    my $type = shift @_ || 'all';

    #--- Start HTML templates ---
    my $TABLE_HTML = qq(
<TABLE cellspacing="0" 
       cellpadding="0" 
       border     ="0" 
       width      ="100%%" 
       class      ="background1">
 <TR>
  <TD colspan="6"><img src="/gfx/EnsMart/blank.gif" height="4" alt=""></TD>
 </TR>
  %s
 <TR>
  <TD colspan="6"><img src="/gfx/EnsMart/blank.gif" height="4" alt=""></TD>
 </TR>
</TABLE>);

    my $VEGA_TABLE_HTML = qq(
<TABLE cellspacing="0" 
       cellpadding="0" 
       border     ="0" 
       width      ="100%%" 
       class      ="background1">
 <TR>
  <TD colspan="8"><img src="/gfx/EnsMart/blank.gif" height="4" alt=""></TD>
 </TR>
  %s
 <TR>
  <TD colspan="8"><img src="/gfx/EnsMart/blank.gif" height="4" alt=""></TD>
 </TR>
</TABLE>);

    my $FORM_HTML = qq(
 <FORM action="%s" method="GET" name="feederform">
 <TR valign="center">
  <TD align=right nowrap>&nbsp; Find %s  </TD>   
  <TD align=right nowrap>&nbsp;      %s  </TD>
  <TD align=right nowrap>&nbsp;      %s  </TD>
  <TD class="barial" nowrap width="100%%">
                         &nbsp;      %s  </TD>
  <TD>                               %s  </TD>
  <TD nowrap>     &nbsp;&nbsp;           </TD>
 </TR>
 </FORM>);

    my $VEGA_FORM_HTML = qq(
    	<form action="%s" method="GET" name="feederform">
    	<tr valign="center">
    	<td align=left nowrap> %s  </td> 
	<td align=left width="100%%"> &nbsp; on %s</td>
    	<td align=right nowrap>&nbsp; Find %s  </td>
    	<td align=right nowrap>&nbsp; %s  </td>
    	<td class="barial" nowrap>  &nbsp; %s  </td>
        <td align=right> &nbsp; </td> 
    	<td align=left nowrap>  %s  </td>
        <td align=right> &nbsp; </td> 
    	</form>);

    my $NO_FORM_HTML = qq(
 <TR valign="center">
  <TD align=right nowrap>&nbsp;</TD>   
  <TD align=right nowrap>&nbsp;</TD>
  <TD align=right nowrap>&nbsp;</TD>
  <TD class="barial" nowrap width="100%%">&nbsp;</TD>
  <TD>                               %s  </TD>
  <TD nowrap>     &nbsp;&nbsp;           </TD>
 </TR>);

    my $SELECT_FIND_HTML = qq(
   <SELECT name="type" class="white" onchange="init()">%s
   </SELECT>);

    my $OPTION_FIND_HTML = qq(
    <OPTION value="%s" %s>%s</OPTION>);

    my $INPUT_FIND_HTML = qq(
   <INPUT name="q" size="24" value="%s">);

    my $SUBMIT_FIND_HTML = qq(
   <INPUT type="image" value="Search" border="0"
          src ="/gfx/lookup.gif">);

    my $SEARCH_EXAMPLE_LINK = qq(<A HREF="%s">%s</A>);

    my $SEARCH_EXAMPLE_LINKS = qq(<B>[e.g. %s]</B>);

    my $SUBMIT_HELP_HTML = qq(
   <A HREF=
    "javascript:window.open('/%s/helpview?se=1&kw=%s',
                            'helpview',
                            'width=400,height=500,resizable,scrollbars');
                void(0);">
    <IMG border="0" align="right" alt="Click for help on %s" 
         src="/gfx/helpview/help.gif">
   </a>);

    #--- End HTML templates ---

    #--- Start Content gathering ---

    #my $PREFIX          = $SPECIES_DEFS->ENSEMBL_PREFIX;
    my $ENSEMBL_SEARCH
        = '';    #$SPECIES_DEFS->get_config('Multi','ENSEMBL_SEARCH');

    my $ENSEMBL_SCRIPT  = lc($ENV{'ENSEMBL_SCRIPT'});

    #put script name in a format for vega header
    my $script_name = $ENSEMBL_SCRIPT;
    $script_name =~ s/.*\/(.*?)$/$1/;
    $script_name = ucfirst($script_name);
    $script_name =~ s/view$/View/;

    my $ENSEMBL_SPECIES = $ENV{ENSEMBL_SPECIES};
    my $help_link       = "";

    my @species = ();    #$SPECIES_DEFS->valid_species( $ENSEMBL_SPECIES );
         #if( ! @species ){ @species = $SPECIES_DEFS->valid_species() }
    if (!@species) { return '' }

    my $form_action = "/${ENSEMBL_SPECIES}/${ENSEMBL_SEARCH}";
    $help_link = "/${ENSEMBL_SPECIES}/helpview?" . "se=1&kw=$ENSEMBL_SCRIPT";

    # Select options
    # Tests for applicability of select options

    my %valid_options = ();
    foreach my $sp (@species) {

        my $idx_ref
            = undef;   #$SPECIES_DEFS->get_config( $sp, 'ENSEMBL_SEARCH_IDXS' );
        ref($idx_ref) eq 'ARRAY' || next;

        # Define index tests
        my $db_ref = undef;    #$SPECIES_DEFS->get_config( $sp, 'databases' );
        ref($db_ref) eq 'HASH' || next;
        my %option_tests = (
            'SNP' => $db_ref->{ENSEMBL_SNP} || $db_ref->{ENSEMBL_GLOVAR},
            'EST' => $db_ref->{ENSEMBL_EST},
            'Disease' => $db_ref->{ENSEMBL_DISEASE}
        );

        # Apply tests and register indexes
        foreach my $idx (@$idx_ref) {
            if (exists($option_tests{$idx})) {
                next if !$option_tests{$idx};
            }
            $valid_options{$idx} = uc($idx) eq uc($type) ? "selected" : undef();
        }
    }

    # Get search links from <species>.ini
    my $SEARCH_LINKS = {};    #$SPECIES_DEFS->SEARCH_LINKS || {};
    my $srch_links;
    my %srch_order;

    # Do we have links for this script? default?
    my %which = map { /^(.+)\d+/i ? (lc($1), 1) : () }
        keys %$SEARCH_LINKS;

    my $use_this;
    if    ($which{$ENSEMBL_SCRIPT}) { $use_this = $ENSEMBL_SCRIPT }
    elsif ($which{'default'})       { $use_this = 'default' }

    foreach my $param (keys %$SEARCH_LINKS) {
        if (my ($name, $type) = ($param =~ /^($use_this\d+)_(TEXT|URL)/i)) {
            my ($num) = ($name =~ /(\d+)/o);
            $srch_links->{$name}->{$type} = $SEARCH_LINKS->{$param};
            $srch_order{$num} = $name;
        }
    }
    my @srch_entries
        = map { $srch_order{$_} } sort { $a <=> $b } keys(%srch_order);

    #--- End Content gathering ---

    #--- Start HTML generation ---

    #Find select box
    my $option_find_html = sprintf($OPTION_FIND_HTML, 'All', '', 'All');
    foreach my $option (sort { $a cmp $b } keys(%valid_options)) {

        $option_find_html .= sprintf($OPTION_FIND_HTML,
            $option, $valid_options{$option} || '', $option);
    }
    my $select_find_html = sprintf($SELECT_FIND_HTML, $option_find_html);

    #Find example links
    my @examples = ();
    my $example_links;
    unless ($ENV{ENSEMBL_SCRIPT} eq 'sitemap') {
    foreach (@srch_entries) {
        push(
            @examples,
            sprintf(
                $SEARCH_EXAMPLE_LINK,
                $srch_links->{$_}->{URL},
                $srch_links->{$_}->{TEXT}
            )
        );
    }
        $example_links = sprintf($SEARCH_EXAMPLE_LINKS, join(', ', @examples));
    }

    #The table row containing the form
    my $form_html;
    if ($ENSEMBL_SEARCH) {

        #      if ($SPECIES_DEFS->SITE_TYPE eq 'Vega') {
        #	  $form_html= sprintf(   $VEGA_FORM_HTML,
        #				 $form_action,
        #				 sprintf($SUBMIT_HELP_HTML,
        #					 $ENSEMBL_SPECIES,
        #					 $ENSEMBL_SCRIPT ),
        #				 $script_name,
        #				 $select_find_html,
        #				 sprintf($INPUT_FIND_HTML, $q),
        #				 $SUBMIT_FIND_HTML,
        #				 @examples ? $example_links : '',
        #			     );
        #      } else {
        $form_html = sprintf($FORM_HTML,
            $form_action,
            $select_find_html,
            sprintf($INPUT_FIND_HTML, $q),
            $SUBMIT_FIND_HTML,
            @examples ? $example_links : '',
            sprintf($SUBMIT_HELP_HTML, $ENSEMBL_SPECIES, $ENSEMBL_SCRIPT));

        #     }
    } else {
        $form_html = sprintf($NO_FORM_HTML,
            sprintf($SUBMIT_HELP_HTML, $ENSEMBL_SPECIES, $ENSEMBL_SCRIPT));
    }
    my $html;

    #  if ($SPECIES_DEFS->SITE_TYPE eq 'Vega') {
    #      $html = sprintf( $VEGA_TABLE_HTML,
    #			  $form_html );
    #  } else {
    $html = sprintf($TABLE_HTML, $form_html);

    #  }
    #--- End HTML generation ---

    return $html;
}

#----------------------------------------------------------------------
sub ensembl_help_link {

    my ($kw, $se) = @_;
    my $LINK = (defined $se ? "se=$se&" : "") . "kw=$kw";
    return
        qq(javascript:X=window.open('/$ENV{'ENSEMBL_SPECIES'}/helpview?$LINK','helpview','height=400,width=500,left=100,screenX=100,top=100,screenY=100,resizable,scrollbars=yes');X.focus();void(0));
}

#----------------------------------------------------------------------
1;
