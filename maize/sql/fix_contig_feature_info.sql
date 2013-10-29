/*
List all contigs with their feature coordinates vs. seq_region coordinates
select att.value,
       f.seq_region_start feature_start, f.seq_region_end feature_end,
       f.seq_region_end - f.seq_region_start + 1 feature_length,
       asm.asm_start, asm.asm_end,
       asm.asm_end - asm.asm_start + 1 asm_length
  from misc_attrib att, misc_feature f, seq_region s, assembly asm
 where f.misc_feature_id = att.misc_feature_id
   and s.name = att.value
   and asm.cmp_seq_region_id = s.seq_region_id
   and att.attrib_type_id = (select attrib_type_id
                               from attrib_type where code = 'name')
   and att.value like 'ctg%'
 order by convert(substr(att.value, 4), decimal) \g
 */
 
/* Set feature coordinates for contigs to be seq_region coordinates */
update misc_feature f set
   f.seq_region_id = (
       select distinct asm.asm_seq_region_id
         from misc_attrib att, seq_region s, assembly asm
        where f.misc_feature_id = att.misc_feature_id
          and s.name = att.value
          and asm.cmp_seq_region_id = s.seq_region_id
          and att.attrib_type_id = 1
          and att.value like 'ctg%'
   ),
   f.seq_region_start = (
       select distinct asm.asm_start
         from misc_attrib att, seq_region s, assembly asm
        where f.misc_feature_id = att.misc_feature_id
          and s.name = att.value
          and asm.cmp_seq_region_id = s.seq_region_id
          and att.attrib_type_id = 1
          and att.value like 'ctg%'
   ),
   f.seq_region_end = (
       select distinct asm.asm_end
         from misc_attrib att, seq_region s, assembly asm
        where f.misc_feature_id = att.misc_feature_id
          and s.name = att.value
          and asm.cmp_seq_region_id = s.seq_region_id
          and att.attrib_type_id = 1
          and att.value like 'ctg%'
   )
   where misc_feature_id in (select misc_feature_id
                               from misc_attrib a
                              where a.attrib_type_id = 1
                                and a.value like 'ctg%');
