delete d2.* from dna_align_feature d1 join dna_align_feature d2 using (seq_region_id, seq_region_start, seq_region_end, seq_region_strand, hit_name, hit_start, hit_end, hit_strand, analysis_id) where d1.dna_align_feature_id < d2.dna_align_feature_id;

