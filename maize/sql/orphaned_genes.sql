delete
  from gene
       left join seq_region using (seq_region_id)
 where seq_region.seq_region_id is null
  