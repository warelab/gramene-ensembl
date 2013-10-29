create temporary table bad_regions as
select * from (
select a.seq_region_id, a.name, v.value as version, a.length,
       min(asm.asm_start) as ctg_start,
       max(asm.asm_end) as ctg_end
  from seq_region a
       left join assembly asm on asm.asm_seq_region_id = a.seq_region_id
       left join seq_region f on f.seq_region_id = asm.cmp_seq_region_id
       left join seq_region_attrib v on v.seq_region_id = a.seq_region_id
                 and v.attrib_type_id = (select attrib_type_id from attrib_type where code = 'acc-version')
   and a.coord_system_id = 1
 group by a.seq_region_id) regions
 where ctg_start != 1 or ctg_end != length;

create temporary table archive as
select seq_region.*
  from seq_region, bad_regions
 where seq_region.name like concat(bad_regions.name, '.%')
   and seq_region.coord_system_id = 1;

-- select seq_region.name, concat(bad_regions.name, '.', bad_regions.version, substr(seq_region.name, 9))
--   from bad_regions, seq_region
--  where seq_region.name like concat(bad_regions.name, '%')
--    and seq_region.name not like '%.%';
-- 
-- select seq_region.name, concat(substr(seq_region.name, 1, 8), substr(seq_region.name from 11))
--   from seq_region, archive
-- where seq_region.name like concat(archive.name, '%');

update seq_region, bad_regions
   set seq_region.name = concat(bad_regions.name, '.', bad_regions.version, substr(seq_region.name, 9))
 where seq_region.name like concat(bad_regions.name, '%')
   and seq_region.name not like '%.%';

update seq_region, archive
   set seq_region.name = concat(substr(seq_region.name, 1, 8), substr(seq_region.name from 11))
 where seq_region.name like concat(archive.name, '%');