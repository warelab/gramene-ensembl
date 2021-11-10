select genome_db.name as SPECIES, seq_member.stable_id as SPLIT_GENE_TRANSCRIPT from split_genes join seq_member using(seq_member_id) join genome_db using(genome_db_id) order by 1,2;
