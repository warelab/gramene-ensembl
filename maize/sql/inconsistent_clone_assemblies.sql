select contig.name,
       clone.name,
       contig.length,
       -- length(dna.sequence) as dna_length,
       assembly.asm_length,
       assembly.asm_start,
       assembly.asm_end
  from (select *,
               asm_end - asm_start + 1 as asm_length
          from assembly) as assembly
       left join seq_region contig
              on contig.seq_region_id = assembly.cmp_seq_region_id
       left join seq_region clone
              on clone.seq_region_id = assembly.asm_seq_region_id
       -- join dna
       --        on dna.seq_region_id = contig.seq_region_id
 where assembly.asm_length != contig.length
   and clone.seq_region_id is not null;
