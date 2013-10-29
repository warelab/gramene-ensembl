package TIGR::Database;

use strict;
use warnings;
use Carp;
use TIGR::Root;
use vars qw(@ISA);

@ISA=qw(TIGR::Root);

my %_attr_data = #	DEFAULT    	ACCESSIBILITY
                  (
			dbh    => [undef, 		1], 
			prefix => ['tigr_', 		1],  #lc for mysql
                    );
# 

#probably not necessary to do anything special
#sub new {
#    my $self =shift;
#    ref($self) or $self=bless {},$self;
#    $self=SUPER::new($self,@_);
#}

sub _required_fields {
    return qw(dbh);
}
sub _ok_fields {
    return qw(dbh prefix);
}
sub _accessible {
    return  $_attr_data{$_[1]}->[1];
}
sub _default {
    print STDERR "database _default($_[1])\n";
    return  $_attr_data{$_[1]}->[0];
}


#Gets all if no args
sub pseudochromosome {
    my($self,@ids)=@_;
    my $prefix=$self->prefix;

    #perhaps should also allow to specify species(one)
    my( $in)=@ids ? " AND ay.chromosome  IN ("
		  .  join( ',', ('?') x scalar(@ids))
		  .")"
	        : "";
    return $self->_objs_from_query('TIGR::PseudoChromosome',
	qq{
	SELECT ay.tigr_database             as tigr_database
	      ,ay.chromosome           as chromosome
	      ,ay.tigr_current_date         as assembly_date
	      ,coords.end3 	       as length
	      ,touched.tigr_date       as seq_date
	      ,rna.rna_genespk	       as rna_genespk
	      ,procod.protein_codingpk as protein_codingpk
	      ,ai.clone_name           as clone_name_a
	      ,hdr.clone_name          as clone_name_h
	      ,hdr.organism            as organism
	      ,hdr.lineage             as lineage
	      ,ai.asmbl_idpcdata       as asmbl_id
	      ,scfd.scaffoldpk         as scaffoldpk
	FROM ${prefix}pseudochromosome pchr
	    ,${prefix}assembly ay
	    ,${prefix}scaffold scfd
	    ,${prefix}header  hdr
	    ,${prefix}gene_list gl
	    ,${prefix}protein_coding procod
	    ,${prefix}rna_genes rna
	    ,${prefix}asmbl_id ai
	    ,${prefix}coordset coords
	    ,${prefix}seq_last_touched touched
	WHERE scfd.pseudochromosomefk=pchr.pseudochromosomepk
	  AND ay.pseudochromosomefk=pchr.pseudochromosomepk
	  AND hdr.assemblyfk=ay.assemblypk
	  AND gl.assemblyfk=ay.assemblypk
	  AND ai.assemblyfk=ay.assemblypk
	  AND coords.assemblyfk=ay.assemblypk
	  AND procod.gene_listfk=gl.gene_listpk
	  AND rna.gene_listfk=gl.gene_listpk
	  AND touched.headerfk=hdr.headerpk
	}.$in
		,@ids);
}

1;
