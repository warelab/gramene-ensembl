package EnsEMBL::Maize::Component::Feedback;

use strict;
use EnsEMBL::Web::Component;

use Maize::BugTracker;
use Readonly;
use CGI;

Readonly my $FORM_NAME => 'feedback_form';
Readonly my $DEBUG_TAG => qq(<div style="color: orange">[DEBUG MODE]</div>);
Readonly my %FORM_ELEMENTS => (
    'refer_to_url' => {
        'type'     => 'URL',
        'required' => 'no',
        'label'    => 'Referring URL:',
        'name'     => 'refer_to_url',
        'order'    => 1,
    },
    'name' => {
        'type'     => 'String',
        'required' => 'yes',
        'label'    => "Your Name:",
        'name'     => 'name',
        'order'    => 2,
    },
    'email' => {
        'type'     => 'String',
        'required' => 'yes',
        'label'    => "Your Email:",
        'name'     => 'email',
        'order'    => 3,
    },
    'organization' => {
        'type'     => 'String',
        'required' => 'no',
        'label'    => "Organization:",
        'name'     => 'organization',
        'order'    => 4,
    },
    'subject' => {
        'type'     => 'String',
        'required' => 'yes',
        'label'    => "Subject:",
        'name'     => 'subject',
        'order'    => 5,
    },
    'category' => {
        'select'   => 'select',
        'type'     => 'DropDown',
        'name'     => 'category',
        'label'    => 'Category:',
        'required' => 'yes',
        'order'    => 6,
    },
    'bodytext' => {
        'type'     => 'Text',
        'required' => 'yes',
        'label'    => 'Questions/Comments',
        'name'     => 'bodytext',
        'order'    => 7,
        'required' => 'no',
    },
    'action' => {
        'type'     => 'Hidden',
        'name'     => 'action',
        'value'    => 'submit',
        'order'    => 8,
        'required' => 'no',
    },
);

sub show_form {
    my $panel    = shift;
    my ($object) = @_;
    my $script   = $object->script;
    my $form     = EnsEMBL::Web::Form->new($FORM_NAME,
        "/@{[$object->species]}/$script", 'get');
    $form->add_element(
        'type' => 'Information',
        'value' =>
            '<div id="feedback_disclaimer"><b>Note:</b> Please provide a <b><em>valid</em></b> email address. Otherwise it would be difficult for us to send you responses. Thanks!</div>'
    );
    my $referer = $object->param('refer_to_url');
    if (!defined($referer)) {
        $referer = $object->referer;
    }

    my $categories = _fetch_categories();

    my $parameters = _initialize_parameters($object);
    $parameters->{'refer_to_url'}->{'value'} = $referer;
    $parameters->{'category'}->{'values'}    = $categories;

    for my $field (
        sort { $parameters->{$a}->{'order'} <=> $parameters->{$b}->{'order'} }
        keys %$parameters
        )
    {
        my %element_arguments = ();
        for my $key (keys %{ $parameters->{$field} }) {
            $element_arguments{$key} = $parameters->{$field}->{$key};
        }
        $form->add_element(%element_arguments);
    }

    if (CGI::param('debug')) {
        $form->add_element(
            'type'  => 'Hidden',
            'name'  => 'debug',
            'value' => 1,
        );
    }

    $form->add_element(
        'type'  => 'Hidden',
        'name'  => 'jumpto_url',
        'value' => $object->referer
    );

    $form->add_button('submit', 'Send your feedback');
    $form->add_button('reset',  'Reset');

    $panel->print($form->render);
    return 1;
}

=pod

=head2 process
    Execute the send operation

=cut

sub process {
    my $panel = shift;
    my ($object) = @_;
    if (CGI::param('debug')) {
        $panel->print($DEBUG_TAG);
    } else {
        $object->send_email;
        $panel->print(<<HTML);
<div id="feedback_success">
    <p>Thank you, <b>@{[$object->param('name')]}</b>.</p>
    <p>Your message was successfully sent to the Maize Helpdesk. We will get back to you shortly.</p>
    <p id="jumpto_link">[<a href="@{[CGI::param('jumpto_url')]}">Return to my last page</a>]</p>
</div>
HTML
        $object->report_bug;
    }
    return 1;
}

sub show_errors {
    my $panel = shift;
    my ($object) = @_;
    if (CGI::param('debug')) {
        $panel->print($DEBUG_TAG);
    }
    if (defined $panel->{_errors}) {
        my @errors = (
            qq(<div id="errors">),
            join("<br/>\n",
                map {"<span class=\"error\">$_</span>"}
                    @{ $panel->{_errors} }),
            qq(</div>)
        );
        $panel->print(@errors);
    }
    return 1;
}

sub is_required {
    my ($param) = @_;
    return 0 unless defined $FORM_ELEMENTS{$param};
    return ($FORM_ELEMENTS{$param}->{'required'} eq 'yes');
}

sub _fetch_categories {
    my $tracker = new Maize::BugTracker;
    my @categories = ({ 'name' => '-', 'value' => q{} });
    for my $category (@{ $tracker->get_categories }) {
        push @categories, +{ 'name' => $category, 'value' => $category };
    }
    return \@categories;
}

sub _initialize_parameters {
    my ($object) = @_;
    my %parameters = ();
    for my $field (keys %FORM_ELEMENTS) {
        $parameters{$field} = {};
        for my $key (keys %{ $FORM_ELEMENTS{$field} }) {
            $parameters{$field}->{$key} = $FORM_ELEMENTS{$field}->{$key};
        }
        $parameters{$field}->{'value'} ||= $object->param($field);
    }
    return \%parameters;
}

1;
