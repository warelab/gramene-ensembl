/*
 * Reconstruct stable IDs by joining the unique portion to the end of the
 * originating clone name. Since the features are mapped to the contig, we
 * first strip the contig portion of the seq_region name and are left with 
 * a versioned accession.
 */

/* GENE STABLE ID */
update gene_stable_id i,
       gene f,
       seq_region s
   set i.stable_id = concat(substring_index(s.name, '-', 1), '_', substring_index(i.stable_id, '_', -1))
 where i.gene_id = f.gene_id
   and f.seq_region_id = s.seq_region_id

/* TRANSCRIPT STABLE ID */
update transcript_stable_id i,
       transcript f,
       seq_region s
   set i.stable_id = concat(substring_index(s.name, '-', 1), '_', substring_index(i.stable_id, '_', -1))
 where i.transcript_id = f.transcript_id
   and f.seq_region_id = s.seq_region_id

/* TRANSLATION STABLE ID
 * Join seq_region via the transcript table
 */
update translation_stable_id i,
       translation f,
       transcript t,
       seq_region s
   set i.stable_id = concat(substring_index(s.name, '-', 1), '_', substring_index(i.stable_id, '_', -1))
 where i.translation_id = f.translation_id
   and f.transcript_id = t.transcript_id
   and t.seq_region_id = s.seq_region_id

/* EXON STABLE ID
 * Join seq_region via the exon_transcript table
 */
update exon_stable_id i,
       exon_transcript f,
       transcript t,
       seq_region s
   set i.stable_id = concat(substring_index(s.name, '-', 1), '_', substring_index(i.stable_id, '_', -1))
 where i.exon_id = f.exon_id
   and f.transcript_id = t.transcript_id
   and t.seq_region_id = s.seq_region_id
