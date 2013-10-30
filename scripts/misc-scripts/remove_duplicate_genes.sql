create temporary table unique_genes as
select analysis_id,
       seq_region_id, seq_region_start, seq_region_end,
       seq_region_strand, min(gene_id) gene_id
  from gene left join analysis using (analysis_id)
 group by 1, 2, 3, 4, 5;

create index unique_genes_idx
       using btree on unique_genes (gene_id);

delete gene
  from gene left join unique_genes using (gene_id)
 where unique_genes.gene_id is null;
 
delete gene_attrib
  from gene_attrib
       left join gene using (gene_id)
 where gene.gene_id is null;

delete gene_stable_id
  from gene_stable_id
       left join gene using (gene_id)
 where gene.gene_id is null;

delete transcript
  from transcript
       left join gene using (gene_id)
 where gene.gene_id is null;

delete transcript_stable_id
  from transcript_stable_id
       left join transcript using (transcript_id)
 where transcript.transcript_id is null;

delete exon_transcript
  from exon_transcript
       left join transcript using (transcript_id)
 where transcript.transcript_id is null;

delete exon
  from exon
       left join exon_transcript using (exon_id)
 where exon_transcript.exon_id is null;

delete exon_stable_id
  from exon_stable_id
       left join exon using (exon_id)
 where exon.exon_id is null;

delete translation
  from translation
       left join transcript using (transcript_id)
 where transcript.transcript_id is null;

delete translation_stable_id
  from translation_stable_id
       left join translation using (translation_id)
 where translation.translation_id is null;

delete protein_feature
  from protein_feature
       left join translation using (translation_id)
 where translation.translation_id is null;