package TIGR::Root;

use strict;
use warnings;
use vars qw( $AUTOLOAD );
use Carp;

use constant PUB_LOCUS_PREFIX => '"LOC_"';

my %_attr_data = #	DEFAULT    	ACCESSIBILITY
                  (
			dbh => [undef, 		1], 
                    );
my %_sth_cache;
# 
sub new {
    my ($caller, %arg)=@_;
    my $caller_is_obj=ref($caller);
    my $class = $caller_is_obj || $caller;
    my $self = bless {}, $class;
#    print STDERR "$class\n";
    for my $attr ( $self->_ok_fields ) {
#	print STDERR " $attr ";
	if(exists $arg{$attr}) {
	    $self->{$attr}=$arg{$attr};
	    delete $arg{$attr};
#	    print STDERR "given";
	} elsif($caller_is_obj && exists $caller->{$attr}) {
	    $self->{$attr}=$caller->{$attr};
#	    print STDERR "copied";
	} elsif ( defined $self->_default($attr)) {
	    $self->{$attr}=$self->_default($attr);
#	    print STDERR "defaulted";
	}
#	print STDERR "\n";
    }
    croak ("new $class invalid attributes ".join(',',keys %arg)) if %arg;
    my @missing= grep { ! exists $self->{$_} }  $self->_required_fields;
    croak ("new $class missing ".join(',',@missing)) if @missing;
    return $self;
}

#will override these three always:
sub _required_fields {
    return qw(dbh);
}
sub _ok_fields {
    return qw(dbh);
}
sub _accessible {
    return  $_attr_data{$_[1]}->[1];
}
sub _default {
    print STDERR "database _default($_[1])\n";
    return  $_attr_data{$_[1]}->[0];
}

#Assuming class name=table name
sub _primary_key {
    my($caller)=@_;
    my $class=ref($caller) || $caller;
    $class=~s/.*://;
    return lc($class."PK");
}
sub _foreign_key {
    my($caller)=@_;
    my $class=ref($caller) || $caller;
    $class=~s/.*://;
    return lc($class."FK");
}

sub _objs_from_query {
    my($self,$newclass,$query,@bindvalues)=@_;

    my $dbh=$self->dbh;

    eval "require $newclass" or croak("require:$@");

    my $sth= ( $_sth_cache{$dbh}{$query}   #if we use many db connections this could accumulate unused statement handles and be a problem
	       ||= $dbh->prepare($query) )
	   or croak ("prepare $query:".$dbh->errstr);
    $sth->execute(@bindvalues) or croak("execute $query (@bindvalues):".$dbh->errstr);
    my @results;
    while(my $row= $sth->fetchrow_hashref('NAME_lc')) {
	push @results,$newclass->new(dbh=>$dbh,prefix=>$self->prefix,parent=>$self,%$row);
    }
    $sth->errstr and croak("fetch $query (@bindvalues):".$sth->errstr);
    $sth->finish;

    return \@results;
}

sub full_description {
    my($self)=@_;
    my @descs;
    push @descs,$self->com_name if $self->_accessible('com_name');
    push @descs,$self->pub_comment if $self->_accessible('pub_comment');
    push @descs,map { $_->{description} } @{$self->xref} if $self->can('xref');
    @descs=sort { length($b) <=> length($a) || $b cmp $a } 
		grep { defined $_ } @descs;
    my $desc = shift @descs;
    for my $d (@descs) {
	last unless $d;
	$desc.="; $d" if index( uc($desc),uc($d))<0; #i.e. append this if don't already have it
    }
    return $desc;
}

sub DESTROY {}  # so it's not handled by AUTOLOAD

sub AUTOLOAD {
    no strict "refs";
    my ($self,$newval)=@_;

    $AUTOLOAD =~ /.*::(\w+)/;
    my $attr=$1;

    if ($self->_accessible($attr)) {

	*{$AUTOLOAD} = sub {
	    if (defined $_[1]) { $_[0]->{$attr} = $_[1] }
	    return $_[0]->{$attr};
	};    ### end of created subroutine

###  this is called first time only
	if (defined $newval) {
	    $self->{$attr} = $newval;
	}
	return $self->{$attr};

    }
    croak "No such method: $AUTOLOAD";

}

1;
