package EnsEMBL::Maize::Document::HTML::Content;
use strict;
use EnsEMBL::Web::Document::HTML::Content;

@EnsEMBL::Maize::Document::HTML::Content::ISA
    = qw(EnsEMBL::Web::Document::HTML::Content);

sub render {
    my $self = shift;
    $self->print(<<HTML);
<div id="container">
    <div id="center" class="column">
HTML
    $self->SUPER::render();
    $self->print(<<HTML);
    </div><!-- center -->
</div><!-- container -->
HTML
    return 1;
}

1;

