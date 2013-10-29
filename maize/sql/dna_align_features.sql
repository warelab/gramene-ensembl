select -- dna_align_feature.*
       count(distinct dna_align_feature.hit_name)
  from seq_region a
       left join assembly asm on asm.asm_seq_region_id = a.seq_region_id
       left join seq_region f on f.seq_region_id = asm.cmp_seq_region_id
       left join seq_region_attrib v on v.seq_region_id = a.seq_region_id
                 and v.attrib_type_id = (select attrib_type_id from attrib_type where code = 'acc-version')
       left join dna_align_feature on dna_align_feature.seq_region_id = f.seq_region_id
       left join analysis on analysis.analysis_id = dna_align_feature.analysis_id
 where a.name = 'AC211704.1'
   and a.coord_system_id = 1
   and analysis.logic_name = 'Maize_est'
 -- group by dna_align_feature.analysis_id;
