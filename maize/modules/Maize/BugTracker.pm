package Maize::BugTracker;

=head1 NAME

Maize::BugTracker

=head1 SYNOPSIS

  use Maize::BugTracker

  my $bug_tracker = new Maize::BugTracker;
  $bug_tracker->set_project('Public Maize Feedback');
  my @categories = @{ $bug_tracker->get_categories };
  $bug_tracker->complain({
      'summary'     => 'Fix this module',
      'description' => 'This module sucks',
      'category'    => 'Modules',
  });
  
=head1 DESCRIPTION

A generic interface for programmatic submission of bugs/issues into the bug-tracker (specifically, Mantis). Uses SOAP to contact the MantisConnect
API on the server to create the bug.

A project is required to be associated with the module, though a default project

=head1 SEE ALSO

    <a href="http://www.mantisbt.org">Mantis Bug Tracker</a>
    <a href="http://futureware.biz/mantisconnect">MantisConnect</a>
    SOAP::Lite

=cut

use strict;
use warnings;

use SOAP::Lite;

use Readonly;

Readonly my $SOAP_PROXY     => 'http://ware.cshl.edu/bugs/mc/mantisconnect.php';
Readonly my $SOAP_NAMESPACE => 'http://futureware.biz/mantisconnect';

Readonly my $MANTIS_USERNAME => 'feedback';
Readonly my $MANTIS_PASSWORD => 'f33db@ck3';

Readonly my $DEFAULT_PROJECT => 'Public Maize Feedback';

=pod

=head2 new
    Constructor

=cut

sub new {
    my $class    = shift;
    my ($params) = @_;
    my $self     = bless {}, $class;
    $self->set_project($params->{'project'} || $DEFAULT_PROJECT);
    return $self;
}

=pod

=head2 init_soap
    Initializes the SOAP object

=cut

sub soap {
    my $self = shift;
    if (!defined $self->{'soap'}) {
        $self->{'soap'} = SOAP::Lite->uri($SOAP_NAMESPACE)->proxy($SOAP_PROXY);
    }
    return $self->{'soap'};
}

=pod

=head2 _invoke_auth
    Wraps around the method to ensure it could be called

=cut

sub _invoke_auth {
    my $self = shift;
    my ($method_name, @arguments) = @_;
    my $soap = $self->soap;
    my $call
        = $soap->$method_name($MANTIS_USERNAME, $MANTIS_PASSWORD, @arguments);
    if ($call->fault) {
        die join('', "Error ", $call->faultcode, ": ", $call->faultstring);
    }
    return $call->result;
}

=pod

=head2 set_project
    Finds an accessible Mantis project and sets it. Dies if project cannot be
    found.

=cut

sub set_project {
    my $self = shift;
    my ($project_name) = @_;

    my @projects = @{ $self->_invoke_auth('mc_projects_get_user_accessible') };
    for my $project (@projects) {
        if ($project->{'name'} eq $project_name) {
            $self->{'_project'} = $project;
        }
    }
    if (!defined $self->{'_project'}) {
        die "Unable to find BugTracker project '$project_name'";
    }
}

=pod

=head2 get_project
    Returns the project

=cut

sub get_project {
    return $_[0]->{'_project'};
}

=pod

=head2 complain
    Creates a new issue

=cut

sub complain {
    my $self = shift;
    my ($params) = @_;

    my @issue_data = (
        [ 'summary',     $params->{'subject'},     'string' ],
        [ 'description', $params->{'description'}, 'string' ],
        [ 'category',    $params->{'category'},    'string' ],
        [ 'project',     $self->get_project,       'ProjectData' ]
    );

    # TODO: This needs to be a lot simpler
    my @array = map { SOAP::Data->name($_->[0] => $_->[1])->type($_->[2]) }
        @issue_data;

    my $issue = SOAP::Data->name(
        'IssueData' => \SOAP::Data->value(
            SOAP::Data->value(@array)->type('IssueData')
        )
    )->type('IssueData');

    my $result = $self->_invoke_auth('mc_issue_add', $issue);
    warn("Issue #$result submitted to Mantis");
}

=pod

=head2 get_categories
    Returns a list of categories for the project

=cut

sub get_categories {
    my $self = shift;
    return $self->_invoke_auth('mc_project_get_categories',
        $self->get_project->{'id'});
}

1;

=pod

=head1 AUTHOR

Shiran Pasternak E<lt>shiran@cshl.eduE<gt>
Apurva Narechania E<lt>apurva@cshl.eduE<gt>

=cut

=head1 COPYRIGHT

Copyright (c) 2006 Cold Spring Harbor Laboratory

This library is free software;  you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut