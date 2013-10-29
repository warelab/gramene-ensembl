package TIGR::Exon;

use strict;
use warnings;
use Carp;
use TIGR::Root;
use vars qw(@ISA);

@ISA=qw(TIGR::Root);

my %_attr_data = #	DEFAULT    	ACCESSIBILITY
                  (
			dbh => [undef, 		1], 
			prefix => ['', 		1], 
			parent => ['', 		1], 
			exonpk => [ undef,	1],
			feat_name => [ undef,	1],
			exon_date => [ undef,	1],
			e_start => [ undef,	1],
			e_end => [ undef,	1],
			coding_start => [ undef,	1],
			coding_end => [ undef,	1],
                    );
# 

sub _required_fields {
    return qw(dbh);
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
