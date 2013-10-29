select acc_att.value as accession,
       nam_att.value as clone_name,
       sta_att.value as status,
       chr.name as chr, mf.seq_region_start as chr_start, mf.seq_region_end as chr_end,
       fpc_att.value as contig,
       mf.seq_region_start - asm.asm_start + 1 as contig_start,
       mf.seq_region_end   - asm.asm_start + 1 as contig_end   
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
 order by chr, chr_start;
