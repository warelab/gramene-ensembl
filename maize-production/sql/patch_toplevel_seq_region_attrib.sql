# Determines all top-level seq_regions in the database and flags these
# as toplevel in the seq_region_attrib table

replace into attrib_type (code, name, description)
      values ('toplevel',
              'Top Level',
              'Top Level Non-Redundant Sequence Region');

delete sa
  from seq_region_attrib sa
       left join attrib_type at using (attrib_type_id)
       left join assembly a on sa.seq_region_id = a.cmp_seq_region_id
 where at.code = 'toplevel'
   and a.cmp_seq_region_id is not null;

replace into seq_region_attrib
select distinct(sr.seq_region_id), at.attrib_type_id, 1
  from seq_region sr
       left join assembly a
              on sr.seq_region_id = a.cmp_seq_region_id,
       attrib_type at
 where cmp_seq_region_id is null
   and at.code = 'toplevel';
