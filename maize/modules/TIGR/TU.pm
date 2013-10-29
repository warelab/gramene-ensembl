package TIGR::TU;

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
			tupk => [undef,	1],
			tu_date => [undef,	1],
			feat_name => [undef,	1],
			name_curated => [undef,	1],
			com_name => [undef,	1],
			pub_comment => [undef,	1],
			locus => [undef,	1],
			ec_num => [undef,	1],
			is_pseudogene => [undef,1],
			alt_locus => [undef,	1],
			pub_locus => [undef,	1],
			gene_sym => [undef,	1],
			tigr_comment => [undef,	1],
			chromo_link => [undef,	1],
			evidence_type => [undef,1],
			assign_method => [undef,1],
			assign_acc => [undef,	1],
			strand =>     [ 1,	1],
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

# returns a list of { db => , description => , id => }
sub xref {
    my ($self)=@_;
    unless( exists $self->{xrefs}) {
        my $rval=[];
	push @$rval,{ db => 'TIGR_TU', id => $self->chromo_link } if $self->chromo_link;
	push @$rval,{ db => 'TIGR_FN', id => $self->feat_name } if $self->feat_name;
	if($self->assign_acc) { #EGAD:, GP:, OMNI:
	    my $db;
	    $self->assign_acc =~ /^([A-Z]+):(\w+):(.*)/
	        and $db=$1 eq 'SP' ? 'SWISSPROT' : $1 
		and push @$rval,{ db => $db, id => $2, description => $3  } ;
	    $self->assign_acc =~ /^(([A-Z]+)\d+):(.*)/
	        and push @$rval,{ db => $2, id => $1, description => $3  } ;
	}
	push @$rval,{ db => 'EC', id => $self->ec_num } if $self->ec_num;

	my %locus;
	push @$rval, map { { db => 'TIGR_LOCUS', id => $_ } } 
		     grep {$_ && !$locus{$_}++} 
		     ( $self->locus , $self->pub_locus , $self->alt_locus );
	$self->{xrefs}=$rval;
    }
    return $self->{xrefs};
}

# Returns a list of gene ontology xrefs
sub go_xref {
  my $self = shift;
  unless( exists $self->{go_xrefs}) {
    my $q = qq{
SELECT  go.assignment,
        go.go_term,
        ev.evidence
FROM    tigr_gene_info info, 
        tigr_gene_ontology ont, 
        tigr_go_id go 
        LEFT JOIN tigr_go_evidence ev ON go.go_idpk = ev.go_idfk
WHERE   info.gene_infopk    = ont.gene_infofk 
AND     ont.gene_ontologypk = go.gene_ontologyfk
AND     info.tufk = ? };
    
    # Should probably handle in TIGR::Root, but would need Xref class
    my $sth = $self->dbh->prepare( $q );
    my $rv = $sth->execute( $self->tupk ) || die( $sth->errstr );
    $self->{go_xrefs} = [];
    foreach my $res( @{$sth->fetchall_arrayref} ){
      push @{$self->{go_xrefs}}, {
        id          => $res->[0],
        description => $res->[1],
        exidence    => $res->[2] || 'IEA',
        db          => 'GO',
      };
    }
  }
  return $self->{go_xrefs} || [];
}

#models, which correspond to transcripts
#Gets all
# NB sometimes >1 chromo_link1 per model, so don't get here
sub model {	
    my($self)=@_;
    my $prefix=$self->prefix;

    return $self->_objs_from_query('TIGR::Model',
	qq{
	SELECT m.modelpk          as modelpk
	      ,m.tigr_date        as model_date
	      ,m.curated          as curated
	      ,CONCAT(}.TIGR::Root::PUB_LOCUS_PREFIX.
	      qq{,m.pub_locus)        as pub_locus
	      ,m.tigr_comment     as tigr_comment
	      ,m.feat_name        as feat_name
	      ,sign(coords.end3-coords.end5)  as strand
	FROM ${prefix}model                   m
	    ,${prefix}coordset             coords
	WHERE m.tufk=?
	  AND coords.modelfk=m.modelpk
	}
	,$self->tupk
    );
}

1;
