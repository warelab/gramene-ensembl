select genome_db.name, seq_member.stable_id from split_genes join seq_member using(seq_member_id) join genome_db using(genome_db_id);
