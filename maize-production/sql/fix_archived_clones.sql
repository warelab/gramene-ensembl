create temporary table ctgs as
select f.seq_region_id as ctg_seq_region_id,
       f.name as ctg_name,
       a.seq_region_id as bad_seq_region_id,
       a.name as bad_name,
       c.seq_region_id as good_seq_region_id,
       c.name as good_name
  from (select * from seq_region where coord_system_id = 1) a
       left join assembly asm on asm.asm_seq_region_id = a.seq_region_id
       left join seq_region f on f.seq_region_id = asm.cmp_seq_region_id
       left join seq_region_attrib v on v.seq_region_id = a.seq_region_id
                 and v.attrib_type_id = (select attrib_type_id
                                           from attrib_type
                                          where code = 'acc-version')
       left join seq_region c on c.name = substr(f.name, 1, 8)
 where f.name not like concat(a.name, '%');
update assembly, ctgs
   set assembly.asm_seq_region_id = ctgs.good_seq_region_id
 where assembly.cmp_seq_region_id = ctgs.ctg_seq_region_id;
 

-- Set the version of seq_version_attrib for archived clones to the
-- version in the seq_region name
update seq_region r, seq_region_attrib a
   set a.value = substr(r.name, 10)
 where a.seq_region_id = r.seq_region_id
   and a.attrib_type_id = (select attrib_type_id
                             from attrib_type
                            where code = 'acc-version')
   and length(r.name) >= 10;

-- For current accessions where the version attrib is null, set the version
-- to the max version of an archived clone and add 1
insert into seq_region_attrib (seq_region_id, attrib_type_id, value)
select current.seq_region_id,
       archive_version.attrib_type_id,
       max(archive_version.value) + 1
 from (select * from seq_region where coord_system_id = 1) current
      left join seq_region_attrib current_version
             on current_version.seq_region_id = current.seq_region_id
            and current_version.attrib_type_id =
                (select attrib_type_id from attrib_type
                  where code = 'acc-version')
      left join (select * from seq_region where coord_system_id = 1) archive
             on archive.name like concat(current.name, '%')
            and length(archive.name) >= 10
      left join seq_region_attrib archive_version
             on archive_version.seq_region_id = archive.seq_region_id
            and archive_version.attrib_type_id =
                (select attrib_type_id from attrib_type
                  where code = 'acc-version')
 where current_version.value is null
 group by current.seq_region_id;


select count(*)
  from seq_region where coord_system_id = 1;
select count(*)
  from seq_region left join seq_region_attrib using (seq_region_id)
 where seq_region_attrib.value = 'current';
select substr(name, 1, 8), count(*)
  from seq_region
 where coord_system_id = 1
 group by substr(name, 1, 8)
having count(*) > 2


select *
  from 