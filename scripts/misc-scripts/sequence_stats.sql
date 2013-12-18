-- set echo on
-- set newpage none
-- set pagesize 3000
-- set sqlselect now() + '""';

--  select 'current time:', CURRENT_TIMESTAMP(); 
--  select '1. Clones by site';

--  col site format A15
--  break on site 

--  select site,htg_phase,count(*) from bacpac group by site,htg_phase ;

--  select '2. Clones by phase';
--  select htg_phase,count(*) from bacpac group by htg_phase;

select '1. Sequence size';
select min(name),sum(asm_end-asm_start) from assembly join seq_region
 where asm_seq_region_id = seq_region_id
  group by seq_region_id with rollup;


select  "2. Contig totals";
select coord_system.name,count(*),sum(length),min(length),max(length) 
from seq_region,coord_system
where seq_region.coord_system_id=coord_system.coord_system_id
      group by coord_system.name;


--  select  '5. Contigs by Chromosome';
--  select chrnum,count(*),sum(length) from bacpac group by chrnum;

--  would need to change for schema change, but don't represent contigs now: */
--  select now() + 'How many Clones have a certain number of pieces';
--  select pieces,count(*) from ( select count(*) pieces,clone_id cid from contig group by clone_id) cc_count group by pieces;

select  '3. Total Genes';
select count(*) from gene;
select  '4. Genes by Analysis';
--  col program format A15
select program ,ngene from analysis
,(select analysis_id aid,count(*) ngene from gene group by analysis_id) gene_by_aid
where analysis.analysis_id=aid;



select  '5. Genes by type';
select biotype,count(*) from gene group by biotype;

select '6. Number of Transcripts';
select count(*) from transcript;

select '7. Number of Translations';
select count(*) from translation;
--  select now() + 'Translation ok? -- 2=ok, 3=Ensembl wrong, 4=Genbank wrong';
--  select ok,count(*) from translation_ok group by ok;

select  '8. Number of Exons';
select count(*) from exon;


select  '9. Exon minimum/average/maximum length and count by chromosome and clone';
select coord_system.name
      ,min(seq_region_end-seq_region_start+1)
      ,avg(seq_region_end-seq_region_start+1)
      ,max(seq_region_end-seq_region_start+1) 
      ,count(*)
from exon,seq_region ,coord_system
where exon.seq_region_id=seq_region.seq_region_id 
   and seq_region.coord_system_id= coord_system.coord_system_id
group by seq_region.coord_system_id;

--  select rank,count(*) from exon_transcript group by rank;

select '10. Count of transcripts by number of exons';
select exons,count(*) 
from ( select count(*) as exons,transcript_id 
       from exon_transcript group by transcript_id) as exon_by_script group by exons;

select  '11. Number of external cross references for each type of object';
--  col ensembl_object_type format A20
select ensembl_object_type ,count(*) 
from object_xref group by ensembl_object_type ;

select  '12. Number of external cross references by destination';
select db_name,xcount 
from external_db,(select external_db_id xdb,count(*) xcount 
                  from xref group by external_db_id) xref_by_ext where xdb=external_db.external_db_id;
--  select externaldbid,count(*) from xref group by externaldbid;

select  '13. Prediction Transcripts';
select logic_name,gff_feature,gff_source,program ,npt
from analysis
    ,(select analysis_id aid,count(*) npt 
      from prediction_transcript group by analysis_id) pt_by_aid
where analysis.analysis_id=aid;

select  '14. Feature hits by type, fcount is number of features, hcount is number of hits';
select logic_name,fcount,hcount from analysis,
(select analysis_id aid,count(distinct hit_name) fcount,count(*) hcount 
 from dna_align_feature group by analysis_id) daf_by_aid 
where analysis.analysis_id=aid
order by logic_name;
-- analysis.db is null for Kiran's stuff

select '15. Marker hits by type, mcount=markers, hcount=hits';
--  col program format A15
--  col logic_name format A13
--  col mcount format 999999
--  col hcount format 999999
select mtype,logic_name,program,mcount,hcount 
from analysis,(select analysis_id aid
                     ,count(distinct marker_feature.marker_id) mcount,count(*) hcount,type mtype 
               from marker_feature,marker 
	       where marker_feature.marker_id=marker.marker_id group by analysis_id,type) marker_by_aid 
where analysis.analysis_id=aid;

select '16. Repeat Feature hits by type, count=consensuses, hcount=hits';
select logic_name,program,ccount,hcount 
from analysis,(select analysis_id aid,count(distinct repeat_consensus_id) ccount,count(*) hcount 
               from repeat_feature 
	       group by analysis_id
              ) repeat_by_aid 
where analysis.analysis_id=aid;


