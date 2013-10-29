package EnsEMBL::Maize::Document::HTML::SearchBox;
use strict;
use CGI qw(escapeHTML);
use EnsEMBL::Web::Document::HTML;
use EnsEMBL::Web::SpeciesDefs;

our $SD = EnsEMBL::Web::SpeciesDefs->new();

our @ISA = qw(EnsEMBL::Web::Document::HTML);

sub new { return shift->SUPER::new('links' => [], 'indexes' => []); }

sub add_link {
    push @{ $_[0]->{'links'} },
        sprintf(
        '<a href="/%s/%s">%s</a>',
        $ENV{'ENSEMBL_SPECIES'},
        CGI::escapeHTML($_[1]),
        CGI::escapeHTML($_[2])
        );
}
sub add_index { push @{ $_[0]->{'indexes'} }, $_[1]; }

sub sp_common : lvalue {
    $_[0]{'sp_common'} =~ s/_/ /g;
    $_[0]{'sp_common'};
}

sub render {
    my $self    = shift;
    my $species = $ENV{'ENSEMBL_SPECIES'};
    if ($species eq 'Multi') { $species = ''; }
    my $common = $species ? $self->sp_common : '';
    my $site_section;
    my $ebang
        = qq(<i><span style="color:#3366bb">e</span><span class="red">!</span></i>);
    if ($common) {
        if ($common =~ /\./) {    # looks like a Latin name
            $site_section = "$ebang <i>$common</i>";
        } else {
            $site_section = "$ebang $common";
        }
    } else {
        $site_section = 'all Ensembl';
    }
    my $script = $SD->ENSEMBL_SEARCH;
    $self->print(
        qq(
<div id="search">
  <form action="/@{[$species||'perl']}/$script" method="get" style="font-size: 0.9em"><div>
    <input type="hidden" name="species" value="@{[$species||'all']}" />
    Search $site_section:
    <select name="idx" style="font-size: 0.9em">
      <option value="All">Anything</option>@{[ map {qq(\n      <option value="$_">$_</option>)} @{$self->{'indexes'}} ]}
    </select>
    <input name="q" size="20" value="@{[CGI::param('q')]}" />
    <input type="submit" value="Go" class="red-button" /></div>
    </form>)
    );
    if (@{ $self->{'links'} }) {
        $self->print(qq(\n  <p class="right" style="margin-right:1em">e.g. ),
            join(", ", @{ $self->{'links'} }), '</p>');
    }
    $self->print(
        qq(
</div>)
    );
}

1;

