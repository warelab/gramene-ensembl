select g.name, gmq.gene_member_stable_id, sm.stable_id from gene_member_qc gmq join seq_member sm using(seq_member_id) join genome_db g on gmq.genome_db_id=g.genome_db_id order by 1,2;


#select genome_db.name as SPECIES, seq_member.stable_id as SPLIT_GENE_TRANSCRIPT from split_genes join seq_member using(seq_member_id) join genome_db using(genome_db_id) order by 1,2;
