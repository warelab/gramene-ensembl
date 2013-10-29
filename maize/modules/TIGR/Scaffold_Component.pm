package TIGR::Scaffold_Component;

use strict;
use warnings;
use Carp;
use TIGR::Root;
use vars qw(@ISA);

@ISA=qw(TIGR::Root);

my %_attr_data = #	DEFAULT    	ACCESSIBILITY
                  (
			dbh => [undef, 		1], 
			prefix => ['TIGR_', 	1], 
			parent => ['', 		1], 
			orientation => [1,	1],
			asmbl_left => [1,	1],
			asmbl_right => [undef,	1],
			chr_left => [1,	1],
			chr_right => [undef,	1],
			id => [undef,	1],
			clone => [undef,	1],
                    );
# 

sub _required_fields {
    return qw(dbh prefix parent chr_left chr_right orientation);
}
sub _ok_fields {
    return keys %_attr_data;
}
sub _accessible {
    return  $_attr_data{$_[1]}->[1];
}
sub _default {
    print STDERR "database _default($_[1])\n";
    return  $_attr_data{$_[1]}->[0];
}

1;
