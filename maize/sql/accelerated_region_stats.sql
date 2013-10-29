create temporary table accelerated_bacs as
select acc_att.value as accession
  from misc_feature mf
       left join misc_attrib nam_att
              on mf.misc_feature_id = nam_att.misc_feature_id
             and nam_att.attrib_type_id = (select attrib_type_id from attrib_type where code = 'name')
       left join misc_attrib fpc_att
              on mf.misc_feature_id = fpc_att.misc_feature_id
             and fpc_att.attrib_type_id = (select attrib_type_id from attrib_type where code = 'superctg')
       left join misc_attrib acc_att
              on mf.misc_feature_id = acc_att.misc_feature_id
             and acc_att.attrib_type_id = (select attrib_type_id from attrib_type where code = 'embl_acc')
       left join misc_attrib ext_att
              on mf.misc_feature_id = ext_att.misc_feature_id
             and ext_att.attrib_type_id = (select attrib_type_id from attrib_type where code = 'external')
       left join misc_attrib sta_att
              on mf.misc_feature_id = sta_att.misc_feature_id
             and sta_att.attrib_type_id = (select attrib_type_id from attrib_type where code = 'seqstatus')
       left join seq_region chr
              on chr.seq_region_id = mf.seq_region_id
       left join seq_region fpc
              on fpc.name = fpc_att.value
       left join assembly asm
              on asm.asm_seq_region_id = chr.seq_region_id
             and asm.cmp_seq_region_id = fpc.seq_region_id
 where ext_att.value = 'false'
   and fpc_att.value in ('ctg181', 'ctg182', 'ctg183')
 order by chr.name, mf.seq_region_start;
 
select g.biotype, count(*)
  from accelerated_bacs ab
       left join zea_mays_core_46_bac.seq_region b
              on ab.accession = b.name
       left join zea_mays_core_46_bac.seq_region_attrib ba
              on ba.seq_region_id = b.seq_region_id
             and ba.attrib_type_id = (select attrib_type_id from zea_mays_core_46_bac.attrib_type where code = 'current-version')
       left join zea_mays_core_46_bac.assembly a
              on a.asm_seq_region_id = b.seq_region_id
       left join zea_mays_core_46_bac.seq_region c
              on c.seq_region_id = a.cmp_seq_region_id
       left join zea_mays_core_46_bac.gene g
              on g.seq_region_id = c.seq_region_id
       left join zea_mays_core_46_bac.gene_stable_id s
              on s.gene_id = g.gene_id
 where ba.value = 'current'
   and s.stable_id is not null
 group by g.biotype
