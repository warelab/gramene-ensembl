package EnsEMBL::Maize::Document::HTML::MastHead;

use strict;
use EnsEMBL::Web::Document::HTML::MastHead;

@EnsEMBL::Maize::Document::HTML::MastHead::ISA
    = qw(EnsEMBL::Web::Document::HTML::MastHead);

sub render {
    my $self = shift;
    my $species_text = '<span style="font-size: 1.5em; color:#fff">&nbsp;</span>';
    if ($self->sp_bio && $ENV{'ENSEMBL_SPECIES'} && $self->sub_title) {
        $species_text
            = qq(<span class="viewname serif">@{[$self->sub_title]}</span>);
    }
    my $masthead = qq{
<div id="masthead">
    <a class="invisible" href="/">
        <div id="maize_logo">&nbsp;</div>
    </a>
    <h1>&nbsp;%s</h1>
</div>
    };
    $self->printf($masthead, $species_text);
}

1;

