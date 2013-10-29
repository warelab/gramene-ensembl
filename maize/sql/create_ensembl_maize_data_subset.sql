delete from local_bac_ensembl.meta;
insert into local_bac_ensembl.meta
select * from meta;

delete from local_bac_ensembl.analysis;
insert into local_bac_ensembl.analysis
select * from analysis;

delete from local_bac_ensembl.attrib_type;
insert into local_bac_ensembl.attrib_type
select * from attrib_type;

delete from local_bac_ensembl.meta_coord;
insert into local_bac_ensembl.meta_coord
select * from meta_coord;

delete from local_bac_ensembl.coord_system;
insert into local_bac_ensembl.coord_system
select * from coord_system;

delete from local_bac_ensembl.seq_region;
insert into local_bac_ensembl.seq_region
select * from (
select seq_region.*
  from seq_region
       left join seq_region_attrib using (seq_region_id)
 where seq_region_attrib.value = 'current'
 limit 10) sr;

delete from local_bac_ensembl.assembly;
insert into local_bac_ensembl.assembly
select * from (
select assembly.*
  from local_bac_ensembl.seq_region
       left join assembly on assembly.asm_seq_region_id = local_bac_ensembl.seq_region.seq_region_id
) a;
       
insert into local_bac_ensembl.seq_region
select * from (
select seq_region.*
  from local_bac_ensembl.assembly
       left join seq_region on seq_region.seq_region_id = local_bac_ensembl.assembly.cmp_seq_region_id
) sr;
 
delete from local_bac_ensembl.seq_region_attrib;
insert into local_bac_ensembl.seq_region_attrib
select * from (
 select seq_region_attrib.*
   from local_bac_ensembl.seq_region
        left join seq_region_attrib using (seq_region_id)
) sra;

delete from local_bac_ensembl.dna;
insert into local_bac_ensembl.dna
select * from (
select dna.*
  from local_bac_ensembl.seq_region
       left join dna using (seq_region_id)
 where dna.seq_region_id is not null
) d;

delete from local_bac_ensembl.gene;
insert into local_bac_ensembl.gene
select * from (
select gene.*
  from local_bac_ensembl.seq_region
       left join gene using (seq_region_id)
 where gene.seq_region_id is not null
) gene;

delete from local_bac_ensembl.transcript;
insert into local_bac_ensembl.transcript
select * from (
select transcript.*
  from local_bac_ensembl.seq_region
       left join transcript using (seq_region_id)
 where transcript.seq_region_id is not null
) transcript;

delete from local_bac_ensembl.simple_feature;
insert into local_bac_ensembl.simple_feature
select * from (
select simple_feature.*
  from local_bac_ensembl.seq_region
       left join simple_feature using (seq_region_id)
 where simple_feature.seq_region_id is not null
) simple_feature;

delete from local_bac_ensembl.repeat_feature;
insert into local_bac_ensembl.repeat_feature
select * from (
select repeat_feature.*
  from local_bac_ensembl.seq_region
       left join repeat_feature using (seq_region_id)
 where repeat_feature.seq_region_id is not null
) repeat_feature;

delete from local_bac_ensembl.dna_align_feature;
insert into local_bac_ensembl.dna_align_feature
select * from (
select dna_align_feature.*
  from local_bac_ensembl.seq_region
       left join dna_align_feature using (seq_region_id)
 where dna_align_feature.seq_region_id is not null
) dna_align_feature;

delete from local_bac_ensembl.protein_align_feature;
insert into local_bac_ensembl.protein_align_feature
select * from (
select protein_align_feature.*
  from local_bac_ensembl.seq_region
       left join protein_align_feature using (seq_region_id)
 where protein_align_feature.seq_region_id is not null
) protein_align_feature;

delete from local_bac_ensembl.misc_feature;
insert into local_bac_ensembl.misc_feature
select * from (
select misc_feature.*
  from local_bac_ensembl.seq_region
       left join misc_feature using (seq_region_id)
 where misc_feature.seq_region_id is not null
) misc_feature;

delete from local_bac_ensembl.misc_set;
insert into local_bac_ensembl.misc_set
select * from misc_set;

delete from local_bac_ensembl.misc_feature_misc_set;
insert into local_bac_ensembl.misc_feature_misc_set
select * from (
select misc_feature_misc_set.*
  from local_bac_ensembl.misc_feature
       left join misc_feature_misc_set using (misc_feature_id)
) misc_feature_misc_set;

delete from local_bac_ensembl.exon_transcript;
insert into local_bac_ensembl.exon_transcript
select * from (
select exon_transcript.*
  from local_bac_ensembl.transcript
       left join exon_transcript using (transcript_id)
 where exon_transcript.transcript_id is not null
) exon_transcript;

delete from local_bac_ensembl.exon;
insert into local_bac_ensembl.exon
select * from (
select exon.*
  from local_bac_ensembl.exon_transcript
       left join exon using (exon_id)
 where exon.exon_id is not null
) exon;

delete from local_bac_ensembl.translation;
insert into local_bac_ensembl.translation
select * from (
select translation.*
  from local_bac_ensembl.transcript
       left join translation using (transcript_id)
 where translation.transcript_id is not null
) translation;

delete from local_bac_ensembl.protein_feature;
insert into local_bac_ensembl.protein_feature
select * from (
select protein_feature.*
  from local_bac_ensembl.translation
       left join protein_feature using (translation_id)
 where protein_feature.translation_id is not null
) protein_feature;

delete from local_bac_ensembl.gene_stable_id;
insert into local_bac_ensembl.gene_stable_id
select * from (
select gene_stable_id.*
  from local_bac_ensembl.gene
       left join gene_stable_id using (gene_id)
 where gene_stable_id.gene_id is not null
) gene_stable_id;

delete from local_bac_ensembl.transcript_stable_id;
insert into local_bac_ensembl.transcript_stable_id
select * from (
select transcript_stable_id.*
  from local_bac_ensembl.transcript
       left join transcript_stable_id using (transcript_id)
 where transcript_stable_id.transcript_id is not null
) transcript_stable_id;

delete from local_bac_ensembl.translation_stable_id;
insert into local_bac_ensembl.translation_stable_id
select * from (
select translation_stable_id.*
  from local_bac_ensembl.translation
       left join translation_stable_id using (translation_id)
 where translation_stable_id.translation_id is not null
) translation_stable_id;

delete from local_bac_ensembl.exon_stable_id;
insert into local_bac_ensembl.exon_stable_id
select * from (
select exon_stable_id.*
  from local_bac_ensembl.exon
       left join exon_stable_id using (exon_id)
 where exon_stable_id.exon_id is not null
) exon_stable_id;

-- delete from local_bac_ensembl.$1;
-- insert into local_bac_ensembl.$1
-- select * from (
-- select $1.*
--   from local_bac_ensembl.$2
--        left join $1 using (${2}_id)
--  where $1.${2}_id is not null
-- ) $1;
-- 
-- $0
-- 

-- delete from local_bac_ensembl.$1;
-- insert into local_bac_ensembl.$1
-- select * from $1;

