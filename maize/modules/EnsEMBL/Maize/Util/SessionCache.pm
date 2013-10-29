package EnsEMBL::Maize::Util::SessionCache;

use strict;
use warnings;

use CGI::Cookie;
use SiteDefs;

use constant ENSEMBL_SESSION_ID => "ENSEMBL_SESSION_ID";

sub get_session {
    return +{};
=pod
    my $request     = Apache->request;
    my $session_id  = _fetch_session_id($request);
    my $session_dir = join('/', $SiteDefs::ENSEMBL_TMP_DIR, 'sessions');
    my $lock_dir    = join('/', $session_dir, 'locks');
    my %session;
    my $age = 'old';

    # print STDERR "Got cookie with session id $session_id\n";
    eval {

        # print STDERR "Attempting to fetch session $session_id\n";
        tie(%session,
            'Apache::Session::File',
            $session_id,
            {   Directory     => $session_dir,
                LockDirectory => $lock_dir,
                Transaction   => 1,
            }
        );

        # print STDERR "Got session $session_id\n";
    };
    if ($@) {

        # print STDERR "Failed to fetch session $session_id: $@\n";
        warn("Unable to fetch session, returning empty hash: $@");
        return +{};
    }
    if (!defined $session_id) {
        $session_id = $session{_session_id};
        _save_session_id($request, $session_id);
        $age = 'new';
    }

    # print STDERR "Using $age session (ID=$session_id)\n";
    return \%session;
=cut
}

sub _fetch_session_id {
    my ($request) = @_;
    my %cookies = CGI::Cookie->fetch($request);
    return
        defined $cookies{ENSEMBL_SESSION_ID}
        ? $cookies{ENSEMBL_SESSION_ID}->value
        : undef;
}

sub _save_session_id {
    my ($request, $session_id) = @_;
    my $cookie = new CGI::Cookie(
        $request,
        -name  => ENSEMBL_SESSION_ID,
        -value => $session_id
    );
    $cookie->bake($request);
}

1;
