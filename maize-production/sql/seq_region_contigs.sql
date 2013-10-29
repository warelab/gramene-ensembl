select count(*) from (select substr(name, 1, 8), count(*) from seq_region where coord_system_id = 1 group by substr(name, 1, 8) having count(*) = 2) c

select s.seq_region_id, s.name
  from seq_region s
       left join seq_region_attrib sa on s.seq_region_id = sa.seq_region_id
       left join attrib_type t on t.attrib_type_id = sa.attrib_type_id
       and t.code = 'acc-version'
 where s.coord_system_id = 1
   and sa.value is null
 limit 10;

select * from attrib_type; 

select * from clone where accession_number = 'AC191595';

select *
  from gene_stable_id
  left join gene on gene_stable_id.gene_id = gene.gene_id
  left join protein_align_feature
    on protein_align_feature.seq_region_id = gene.seq_region_id
   and protein_align_feature.seq_region_start >= gene.seq_region_start
   and protein_align_feature.seq_region_end <= gene.seq_region_end                                
 where stable_id = 'AC191595.2_FG004';

-- Clones whose contigs don't line up with the clone ends
create temporary table bad_regions as
select * from (
select a.seq_region_id, a.name, v.value, a.length,
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
select *
  from seq_region, bad_regions
 where seq_region.name like concat(bad_regions.name, '.%')
   and seq_region.coord_system_id = 1;
    
select a.seq_region_id, a.name, v.value, a.length,
       f.seq_region_id, f.name, f.length, asm.asm_start, asm.asm_end
   from bad_regions a
        left join assembly asm on asm.asm_seq_region_id = a.seq_region_id
        left join seq_region f on f.seq_region_id = asm.cmp_seq_region_id
        left join seq_region_attrib v on v.seq_region_id = a.seq_region_id
                  and v.attrib_type_id = (select attrib_type_id from attrib_type where code = 'acc-version');

select a.seq_region_id, a.name, count(*)
 from bad_regions a
      left join assembly asm on asm.asm_seq_region_id = a.seq_region_id
      left join seq_region f on f.seq_region_id = asm.cmp_seq_region_id
      left join seq_region_attrib v on v.seq_region_id = a.seq_region_id
                and v.attrib_type_id = (select attrib_type_id from attrib_type where code = 'acc-version')
group by a.seq_region_id
having count(*) = 1;

   

select a.seq_region_id, a.coord_system_id, a.name, v.value, a.length, f.seq_region_id, f.name, f.length, asm.asm_start, asm.asm_end
  from (select * from seq_region where coord_system_id = 1) a
       left join assembly asm on asm.asm_seq_region_id = a.seq_region_id
       left join seq_region f on f.seq_region_id = asm.cmp_seq_region_id
       left join seq_region_attrib v on v.seq_region_id = a.seq_region_id
                 and v.attrib_type_id = (select attrib_type_id from attrib_type where code = 'acc-version')
 where f.name not like concat(a.name, '%')

   
select * from seq_region where coord_system_id = 1 and name like 'AC183312%'
update assembly set asm_seq_region_id = 106904
 where cmp_seq_region_id in (select seq_region_id from seq_region where name like 'AC177838-Contig%')
     
     
update seq_region set name = 'AC196277.3' where name = 'AC196277';
update seq_region set name = concat('AC196277', substr(name, 10))
 where name like 'AC1962772%';
 
select * from seq_region where name like 'AC199772-%'
 
select *
  from clone left join contig on clone.id = contig.clone_id
 where clone.accession_number = 'AC195593'
 
delete from assembly where asm_seq_region_id = (select seq_region_id from seq_region where name = 'AC196277.3');
    
    
select * from input_id_analysis where input_id like '%AC202427%';


select count(*)
  from seq_region a
       left join assembly asm on asm.asm_seq_region_id = a.seq_region_id
       left join seq_region f on f.seq_region_id = asm.cmp_seq_region_id
       left join seq_region_attrib v on v.seq_region_id = a.seq_region_id
 where a.coord_system_id = 1
   and f.coord_system_id = 2
   and v.value = 'current'
 order by a.name, asm.asm_start;
 
select count(*) from seq_region where coord_system_id = 2
