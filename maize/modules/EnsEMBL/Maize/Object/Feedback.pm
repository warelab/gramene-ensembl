package EnsEMBL::Maize::Object::Feedback;

use strict;
use warnings;

use EnsEMBL::Web::Object;
use Mail::Mailer;

use Maize::BugTracker;

our @ISA = qw(EnsEMBL::Web::Object);

sub send_email {
    my $self = shift;

    my @message = map { sprintf("%-16.16s %s", "$_->[0]:", $_->[1]) } (
        [ 'Date',         $self->now ],
        [ 'Name',         $self->param('name') ],
        [ 'Email',        $self->param('email') ],
        [ 'Organization', $self->param('organization') ],
        [ 'Referrer',     $self->param('refer_to_url') ],
        [ 'Category',     $self->param('category') ],
        [ 'User agent',   $ENV{'HTTP_USER_AGENT'} ],
        [ 'Subject',      $self->param('subject') ]
    );

    push @message, qq(
Message:

@{[$self->param('bodytext')]}
    );
    my $mailer = new Mail::Mailer('sendmail', Server => "localhost");
    my $sitetype = ucfirst(lc($self->species_defs->ENSEMBL_SITETYPE))
        || 'Ensembl';
    $mailer->open(
        {   'To'      => $self->species_defs->ENSEMBL_HELPDESK_EMAIL,
            'Subject' => "$sitetype Feedback: @{[$self->param('subject')]}",
            'Cc'   => "@{[$self->param('name')]} <@{[$self->param('email')]}>",
            'From' =>
                "MaizeSequence.org Helpdesk <@{[$self->species_defs->ENSEMBL_HELPDESK_EMAIL]}>",
        }
    );
    print $mailer join "\n", @message;
    $mailer->close();
    warn("Sent message to @{[$self->species_defs->ENSEMBL_HELPDESK_EMAIL]}");
    $self->problem(
        'redirect',
        sprintf(
            "/%s/%s?action=thank_you;ref=%s",
            $self->species,              $self->script,
            CGI::escape($self->referer), CGI::escape($self->param('kw')),
            ''
        )
    );
    return 1;

}

=pod

=head2 report_bug
    Reports a bug to the Bug Tracker (Mantis)

=cut

sub report_bug {
    my $self   = shift;
    my @fields = (
        [ 'Submitter',    $self->param('name') ],
        [ 'Email',        $self->param('email') ],
        [ 'Organization', $self->param('organization') ],
        [ 'Referrer',     $self->param('refer_to_url') ],
        [ 'User agent',   $ENV{'HTTP_USER_AGENT'} ]
    );
    my $description = join("\n",
        "Automatic issue submission on @{[$self->now]}",
        (map {"$_->[0]: $_->[1]"} @fields),
        q{}, 'Comments:', $self->param('bodytext'));
    my $tracker = new Maize::BugTracker;
    $tracker->complain(
        {   'subject'     => $self->param('subject'),
            'category'    => $self->param('category'),
            'description' => $description,
        }
    );
    return 1;
}

=pod

=head2 now
    Returns the current date

=cut

sub now {
    my $self = shift;
    my @T    = localtime();
    sprintf "%04d-%02d-%02d %02d:%02d:%02d", $T[5] + 1900, $T[4] + 1, $T[3],
        $T[2], $T[1], $T[0];
}

1;
