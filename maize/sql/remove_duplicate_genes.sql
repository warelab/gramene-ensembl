-- TEMPORARY TABLES FOR OBJECT IDS --
create temporary table extra_genes as
select max(gene_id) as gene_id
  from gene
 group by seq_region_id, seq_region_start, seq_region_end, seq_region_strand
 having count(*) = 2;

create temporary table extra_transcripts as
select transcript_id
  from transcript, extra_genes
 where extra_genes.gene_id = transcript.gene_id;

create temporary table extra_translations as
select translation_id
  from translation, extra_transcripts
 where extra_transcripts.transcript_id = translation.transcript_id;

create temporary table extra_exons as
select exon_id
  from exon_transcript, extra_transcripts
 where extra_transcripts.transcript_id = exon_transcript.transcript_id;

-- EXON TABLES --
select 'Table: exon';
delete
  from exon
 using exon, extra_exons
 where exon.exon_id = extra_exons.exon_id;

select 'Table: exon_stable_id';
delete
  from exon_stable_id
 using exon_stable_id, extra_exons
 where exon_stable_id.exon_id = extra_exons.exon_id;

select 'Table: exon_transcript';
delete
  from exon_transcript
 using exon_transcript, extra_exons
 where exon_transcript.exon_id = extra_exons.exon_id;

select 'Table: supporting_feature';
delete
  from supporting_feature
 using supporting_feature, extra_exons
 where supporting_feature.exon_id = extra_exons.exon_id;


-- TRANSLATION TABLES --
select 'Table: protein_feature';
delete
  from protein_feature
 using protein_feature, extra_translations
 where protein_feature.translation_id = extra_translations.translation_id;

select 'Table: translation_attrib';
delete
  from translation_attrib
 using translation_attrib, extra_translations
 where translation_attrib.translation_id = extra_translations.translation_id;

select 'Table: translation_stable_id';
delete
  from translation_stable_id
 using translation_stable_id, extra_translations
 where translation_stable_id.translation_id = extra_translations.translation_id;

select 'Table: translation';
delete
  from translation
 using translation, extra_translations
 where translation.translation_id = extra_translations.translation_id;

-- TRANSCRIPT TABLES --
select 'Table: regulatory_factor_coding';
delete
  from regulatory_factor_coding
 using regulatory_factor_coding, extra_transcripts
 where regulatory_factor_coding.transcript_id = extra_transcripts.transcript_id;

select 'Table: transcript_attrib';
delete
  from transcript_attrib
 using transcript_attrib, extra_transcripts
 where transcript_attrib.transcript_id = extra_transcripts.transcript_id;

select 'Table: transcript_stable_id';
delete
  from transcript_stable_id
 using transcript_stable_id, extra_transcripts
 where transcript_stable_id.transcript_id = extra_transcripts.transcript_id;

select 'Table: transcript_supporting_feature';
delete
  from transcript_supporting_feature
 using transcript_supporting_feature, extra_transcripts
 where transcript_supporting_feature.transcript_id = extra_transcripts.transcript_id;

select 'Table: unconventional_transcript_association';
delete
  from unconventional_transcript_association
 using unconventional_transcript_association, extra_transcripts
 where unconventional_transcript_association.transcript_id = extra_transcripts.transcript_id;

select 'Table: transcript';
delete
  from transcript
 using transcript, extra_transcripts
 where transcript.transcript_id = extra_transcripts.transcript_id;


-- GENE TABLES --
select 'Table: alt_allele';
delete
  from alt_allele
 using alt_allele, extra_genes
 where alt_allele.gene_id = extra_genes.gene_id;

select 'Table: gene_attrib';
delete
  from gene_attrib
 using gene_attrib, extra_genes
 where gene_attrib.gene_id = extra_genes.gene_id;

select 'Table: gene_stable_id';
delete
  from gene_stable_id
 using gene_stable_id, extra_genes
 where gene_stable_id.gene_id = extra_genes.gene_id;

select 'Table: regulatory_factor_coding';
delete
  from regulatory_factor_coding
 using regulatory_factor_coding, extra_genes
 where regulatory_factor_coding.gene_id = extra_genes.gene_id;

select 'Table: unconventional_transcript_association';
delete
  from unconventional_transcript_association
 using unconventional_transcript_association, extra_genes
 where unconventional_transcript_association.gene_id = extra_genes.gene_id;

select 'Table: gene';
delete
  from gene
 using gene, extra_genes
 where gene.gene_id = extra_genes.gene_id;

select 'Table: input_id_analysis';
delete
  from input_id_analysis
 using input_id_analysis left join translation
    on input_id_analysis.input_id = translation.translation_id
 where input_id_analysis.input_id_type = 'TRANSLATIONID'
   and translation.translation_id is null;

