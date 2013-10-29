package TIGR::Model;

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
			modelpk => [ undef,	1],
			model_date => [ undef,	1],
			curated => [ undef,	1],
			pub_locus => [ undef,	1],
			tigr_comment => [ undef,	1],
			feat_name => [ undef,	1],
			strand => [ 1,	1],
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

sub chromo_link {  #returns pointer to array of strings
    my ($self)=@_;
    my $prefix=$self->prefix;
    return $self->{chromo_link} ||=
           $self->dbh->selectcol_arrayref("SELECT chromo_link FROM ${prefix}chromo_link WHERE modelfk=?"
    			,{},$self->modelpk);
}

sub coords {
    my ($self)=@_;
    my $prefix=$self->prefix;
    return $self->{coords} ||=
           $self->dbh->selectrow_arrayref("SELECT end5,end3 FROM ${prefix}coordset WHERE modelfk=?"
    			,{},$self->modelpk);
}


# returns a list of { db => , description => , id => }
sub xref {
    my ($self)=@_;
    unless( exists $self->{xrefs}) {
        my $rval=[];
	push @$rval,map { { db => 'TIGR_MODEL', id => $_ } } @{$self->chromo_link} ;
		# -use for links to TIGR site
	push @$rval,{ db => 'TIGR_FN', id => $self->feat_name } if $self->feat_name;
		# -need this for checking against .pep file if nothing else
	push @$rval, { db => 'TIGR_LOCUS', id => $self->pub_locus } 
	    if $self->pub_locus;
	$self->{xrefs}=$rval;
    }
    return $self->{xrefs};
}

#exons
#Gets all
sub exon {	
    my($self)=@_;
    my $prefix=$self->prefix;

    my $rvalue= $self->_objs_from_query('TIGR::Exon',
	qq{
	SELECT e.exonpk       as exonpk
	      ,e.feat_name    as feat_name
	      ,e.tigr_date    as exon_date
	      ,coords.end5    as e_start
	      ,coords.end3    as e_end
	      ,cc.end5        as coding_start
	      ,cc.end3        as coding_end
	FROM ${prefix}coordset coords
	    ,${prefix}exon     e
	       LEFT OUTER JOIN ${prefix}cds c
	           ON  c.exonfk=e.exonpk
	       LEFT OUTER JOIN ${prefix}coordset cc
		   ON cc.cdsfk=c.cdspk
	WHERE e.modelfk=?
	  AND coords.exonfk=e.exonpk
	ORDER BY coords.end5*?
	}
	,$self->modelpk
	,$self->strand
    );
    unless(scalar(@$rvalue) ) { #problem in Oryza Sativa 1.0 data 
				#Use gene model coordinates
				#as cds and TU as exon
    }
    return $rvalue;
}


1;
