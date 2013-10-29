package EnsEMBL::Maize::Util::FileCache;

use strict;
use warnings;

use File::Spec;
use Storable;
use SiteDefs;

=pod

=head2 new
    Constructor

=cut

sub new {
    my ($class_name, $params) = @_;
    my $self = bless +{}, $class_name;
    $self->{'cache_filename'} = File::Spec->catfile($SiteDefs::ENSEMBL_TMP_DIR,
        $params->{'filename'});
    $self->_init_cache;
    return $self;
}

sub print {
    my $self       = shift;
    my (@messages) = @_;
    print STDERR "CACHE [@{[ \$self->{'cached_data'} ]}]: ",
        join('', @messages), "\n";
}

sub _load_cache {
    my $self     = shift;
    my $contents = undef;
    eval {
        $self->{'cached_data'} = retrieve($self->{'cache_filename'});
    };
    if ($@) {
        $self->{'cached_data'} = +{};
    }
    $self->{'file_loaded'} = 1;
}

sub get_value {
    my $self        = shift;
    my ($cache_key) = @_;
    if (!defined($self->{'cached_data'})) {
        if ($self->{'file_loaded'}) {
            $self->{'cached_data'} = +{};
        } else {
            $self->_load_cache;
        }
    }
    return $self->{'cached_data'}->{$cache_key};
}

sub set_value {
    my $self = shift;
    my ($cache_key, $value) = @_;
    if (!defined $self->{'cached_data'}) {
        $self->_load_cache;
    }
    $self->{'cached_data'}->{$cache_key} = $value;
}

sub save {
    my $self = shift;
    store($self->{'cached_data'}, $self->{'cache_filename'});
}

=pod

=head2 _init_cache
    Initializes the cache

=cut

sub _init_cache {
    my $self = shift;
}

sub DESTROY {
    my $self = shift;
}

1;