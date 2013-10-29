package EnsEMBL::Maize::Document::HTML::HelpLink;

use strict;

our @ISA = qw(EnsEMBL::Web::Document::HTML::HelpLink);

sub new {
    my $self = shift->SUPER::new(@_);
    return $self;
}

sub label  :lvalue { $_[0]{'label'} = 'Ensembl Help'; $_[0]{'label'}; }
sub simple :lvalue { $_[0]{'simple'} = 1; $_[0]{'simple'}; }

1;

