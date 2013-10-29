select biotype, count(*)
  from gene
       left join seq_region f on gene.seq_region_id = f.seq_region_id
       left join assembly asm on asm.cmp_seq_region_id = f.seq_region_id
       left join seq_region a on a.seq_region_id = asm.asm_seq_region_id
       left join seq_region_attrib on seq_region_attrib.seq_region_id = a.seq_region_id
where seq_region_attrib.value = 'current'
 group by biotype
