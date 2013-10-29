select a.name as clone, a.length as clone_len,
       f.name as contig, f.length as contig_len,
       asm.asm_start as start, asm.asm_end as end,
       asm.asm_end - asm.asm_start + 1 as asm_len
  from seq_region a
       left join assembly asm on asm.asm_seq_region_id = a.seq_region_id
       left join seq_region f on f.seq_region_id = asm.cmp_seq_region_id
       left join seq_region_attrib v on v.seq_region_id = a.seq_region_id
                 and v.attrib_type_id = (select attrib_type_id from attrib_type where code = 'acc-version')
 where a.name = 'AC199375.3'
   and a.coord_system_id = 1
 order by a.name, asm.asm_start;

select *
  from seq_region
       left join assembly on assembly.cmp_seq_region_id = seq_region.seq_region_id
 where name like 'AC199375.3%'
  
-- select *
--   from seq_region where coord_system_id = 1 and seq_region_id > (select seq_region_id from seq_region where name = 'AC198199.3')
--  order by seq_region_id limit 10
--  
--  seq_region_id   name    coord_system_id length
--  187506  AC198197.3      1       120261
--  187535  AC198655.3      1       168787
--  187545  AC198654.2      1       243025
--  187560  AC198652.2      1       162184
--  187572  AC198649.2      1       206122
--  187589  AC198605.2      1       147534
--  187601  AC198602.3      1       185381
--  187609  AC198596.3      1       124748
--  187616  AC198594.3      1       126496
--  187630  AC198739.3      1       176952