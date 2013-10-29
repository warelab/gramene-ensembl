package TIGR::PseudoChromosome;

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
			parent => [undef,	1], 
			tigr_database => [ undef, 1],
			chromosome => [ undef, 1],
			# -- chromosome = number (no leading zeros)
			assembly_date => [ undef, 1],
			length => [ undef, 1],
			seq_date => [ undef, 1],
			rna_genespk => [ undef, 1],
			protein_codingpk => [ undef, 1],
			clone_name_a => [ undef, 1],
			clone_name_h => [ undef, 1],
			# -- clone_name like "CHR09v07232003"
			organism => [ undef, 1],
			lineage => [ undef, 1],
			asmbl_id => [ undef, 1],
			scaffoldpk => [ undef, 1],
                    );
#

sub _required_fields {
    return qw(dbh chromosome length scaffoldpk protein_codingpk );
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

#Gets all 
sub scaffold_component {
    my($self)=@_;
    my $prefix=$self->prefix;

    return $self->_objs_from_query('TIGR::Scaffold_Component',
	qq{
	SELECT sc.orientation       as orientation
	      ,sc.asmbl_left_coord  as asmbl_left
	      ,sc.asmbl_right_coord as asmbl_right
	      ,sc.chr_left_coord    as chr_left
	      ,sc.chr_right_coord   as chr_right
	      ,ai.asmbl_idpcdata    as id
	      ,ai.clone_name        as clone
	FROM ${prefix}scaffold_component sc
	    ,${prefix}asmbl_id ai
	WHERE sc.scaffoldfk=?
	  AND ai.scaffold_componentfk=sc.scaffold_componentpk
	}
	,$self->scaffoldpk);
}


#transcriptional units, which correspond to coding genes
#Gets all 
#assumes 
#       exactly one funct_annot_evidence per gene_info (dtd says *)
#       at most one assign_acc per funct_annot_evidence (dtd says *)
#       exactly one chromo_link per tu (dtd says *)
# dtd does guarantee exacty one gene_info per tu
	      ## removed by kiran to account change in schema
	      #,gi.alt_locus        alt_locus
	      #,gi.ec_num           ec_num
	      #,gi.gene_sym         gene_sym
sub tu {	
    my($self)=@_;
    my $prefix=$self->prefix;

    return $self->_objs_from_query('TIGR::TU',
	qq{
	SELECT tu.tupk             as tupk
	      ,tu.tigr_date        as tu_date
	      ,tu.feat_name        as feat_name
	      ,cn.curated          as name_curated
	      ,cn.com_namepcdata   as com_name
	      ,gi.pub_comment      as pub_comment
	      ,gi.locus            as locus
	      ,gi.is_pseudogene    as is_pseudogene
	      ,CONCAT(}.TIGR::Root::PUB_LOCUS_PREFIX.
	      qq{,gi.pub_locus)        as pub_locus
	      ,gi.tigr_comment     as tigr_comment
	      ,cl.chromo_link      as chromo_link
	      ,fae.type            as evidence_type
	      ,aa.assign_method    as assign_method
	      ,aa.assign_accpcdata as assign_acc
	      ,sign(coords.end3-coords.end5)  as strand
	FROM ${prefix}tu                   tu
	    ,${prefix}gene_info            gi
	    ,${prefix}chromo_link1         cl
	    ,${prefix}com_name             cn
	    ,${prefix}coordset             coords
	    ,${prefix}funct_annot_evidence fae
	  LEFT OUTER JOIN ${prefix}assign_acc aa 
	      ON fae.funct_annot_evidencepk=aa.funct_annot_evidencefk
	WHERE tu.protein_codingfk=?
	  AND cl.tufk=tu.tupk
	  AND gi.tufk=tu.tupk
	  AND cn.gene_infofk=gi.gene_infopk
	  AND fae.gene_infofk=gi.gene_infopk
	  AND coords.tufk=tu.tupk
	}
	,$self->protein_codingpk
    );
}

1;
