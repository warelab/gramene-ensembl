package EnsEMBL::Maize::Configuration::Feedback;

use strict;
use EnsEMBL::Web::Configuration;
use Email::Valid;
use Readonly;

use EnsEMBL::Maize::Component::Feedback;

our @ISA = qw( EnsEMBL::Web::Configuration );

Readonly my $PAGE_TITLE => 'Site Feedback';

sub configure_feedback {
    my $self       = shift;
    my $object     = $self->{object};
    my $print_form = 1;
    my $panel      = new EnsEMBL::Web::Document::Panel(
        'object'  => $self->{object}
    );
    if ($object->param('action') eq 'submit') {
        my @errors = $self->validate($object);
        if (scalar @errors > 0) {
            $panel->{_errors} = \@errors;
            $print_form = 1;
        } else {
            $panel->add_components(
                qw(process_form EnsEMBL::Maize::Component::Feedback::process)
            );
            $print_form = 0;
        }
    }
    if ($print_form) {
        $panel->add_components(
            qw{
                show_errors EnsEMBL::Maize::Component::Feedback::show_errors
                show_form   EnsEMBL::Maize::Component::Feedback::show_form
                }
        );
    }
    $self->{page}->content->add_panel($panel);
    $self->{page}->title->set($PAGE_TITLE);
}

=pod

=head2 context_menu
    Draw the context menu

=cut

sub context_menu {

    # For now, draw nothing
}

=pod

=head2 validate
    Validate parameters

=cut

sub validate {
    my $self            = shift;
    my ($object)        = @_;
    my @required_params = qw(name email subject category);
    my @errors          = ();
    for my $param ($object->param()) {
        if (EnsEMBL::Maize::Component::Feedback::is_required($param)) {
            my $value = $object->param($param);
            if (!$value || $value eq q{}) {
                push @errors, "Required parameter: $param";
            }
        }
    }
    for my $param (@required_params) {
    }
    my $user_email = $object->param('email');
    if ($user_email ne ''
        && !Email::Valid->address(
            -address => $user_email,
            -mxcheck => 1
        )
        )
    {
        push @errors, "Invalid email: $user_email";
    }
    return @errors;
}

1;
