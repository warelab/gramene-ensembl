create temporary table unique_pdt as
select analysis_id,
       seq_region_id, seq_region_start, seq_region_end,
       seq_region_strand, min(prediction_transcript_id) prediction_transcript_id 
  from prediction_transcript left join analysis using (analysis_id)
 group by 1, 2, 3, 4, 5;

create index unique_pdt_idx
       using btree on unique_pdt (prediction_transcript_id);

delete prediction_transcript 
  from prediction_transcript left join unique_pdt using (prediction_transcript_id)
 where unique_pdt.prediction_transcript_id is null;
 

delete prediction_exon
  from prediction_exon 
       left join prediction_transcript using (prediction_transcript_id)
 where prediction_transcript.prediction_transcript_id is null;

