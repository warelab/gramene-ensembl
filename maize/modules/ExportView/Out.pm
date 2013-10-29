###############################################################################
#
#   Name:	    ExportView::Out;
#
#   Description:    Tie Array to this so 'push' outputs with newline
#
#   History:	    2004-03-23	scs - original version
#
###############################################################################
package ExportView::Out;

use strict;
use Carp;
use Symbol;

sub TIEARRAY {
    my($class,$fh)=@_;
    unless ($fh) {
	#open $fh, ">&STDOUT" or croak("can't dup STDOUT");
	$fh=select;
    }
    my $self=\$fh;
    return bless $self,$class;
}

sub PUSH {
    my $self=shift;
    my $oldfh=select($$self);
    for my $newval (@_) {
	print $newval,"\n";
    }
    select($oldfh);
}

sub FETCHSIZE {
    return 0;
}

1;
