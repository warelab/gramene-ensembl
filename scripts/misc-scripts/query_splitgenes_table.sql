select genome_db.name, gene_member.stable_id, sm.stable_id from split_genes sg join seq_member sm using (seq_member_id) join gene_member using(gene_member_id) join genome_db on(gene_member.genome_db_id=genome_db.genome_db_id);

