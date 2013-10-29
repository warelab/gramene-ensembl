alter table dna_align_feature add constraint unique key (seq_region_id, seq_region_start, seq_region_end, seq_region_strand, hit_name, hit_start, hit_end, hit_strand, analysis_id);

