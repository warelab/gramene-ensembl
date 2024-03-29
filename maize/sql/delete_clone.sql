create temporary table delete_this_clone as
 select seq_region_id from seq_region
 where name like '%.0';

create temporary table orphaned_contigs as
select contig.seq_region_id,
       contig.name
  from assembly
       left join seq_region contig
              on contig.seq_region_id = assembly.cmp_seq_region_id
       left join seq_region clone
              on clone.seq_region_id = assembly.asm_seq_region_id
 where clone.seq_region_id is null;

select count(*) as dna_count
  from dna join orphaned_contigs using (seq_region_id);


delete a, asm, f, v, dna,
       dna_align_feature,
       protein_align_feature,
       protein_feature,
       repeat_feature,
       simple_feature,
       misc_feature,
       misc_feature_misc_set,
       gene, gene_stable_id,
       transcript, transcript_stable_id,
       exon, exon_stable_id,
       exon_transcript,
       translation, translation_stable_id
from delete_this_clone dtc
    left join seq_region a
           on a.seq_region_id = dtc.seq_region_id
    left join assembly asm
           on asm.asm_seq_region_id = a.seq_region_id
    left join seq_region f
           on f.seq_region_id = asm.cmp_seq_region_id
    left join seq_region_attrib v
           on v.seq_region_id = a.seq_region_id
    left join dna
           on dna.seq_region_id = f.seq_region_id
    left join gene
           on gene.seq_region_id = f.seq_region_id
    left join gene_stable_id
           on gene_stable_id.gene_id = gene.gene_id
    left join transcript
           on transcript.seq_region_id = f.seq_region_id
    left join transcript_stable_id
           on transcript_stable_id.transcript_id = transcript.transcript_id
    left join exon_transcript
           on exon_transcript.transcript_id = transcript.transcript_id
    left join exon
           on exon.exon_id = exon_transcript.exon_id
    left join exon_stable_id
           on exon_stable_id.exon_id = exon.exon_id
    left join translation
           on translation.transcript_id = transcript.transcript_id
    left join translation_stable_id
           on translation_stable_id.translation_id = translation.translation_id
    left join dna_align_feature
           on dna_align_feature.seq_region_id = f.seq_region_id
    left join protein_align_feature
           on protein_align_feature.seq_region_id = f.seq_region_id
    left join protein_feature
           on protein_feature.translation_id = translation.translation_id
    left join repeat_feature
           on repeat_feature.seq_region_id = f.seq_region_id
    left join simple_feature
           on simple_feature.seq_region_id = f.seq_region_id
    left join misc_feature
           on misc_feature.seq_region_id = f.seq_region_id
    left join misc_feature_misc_set
           on misc_feature_misc_set.misc_feature_id = misc_feature.misc_feature_id
