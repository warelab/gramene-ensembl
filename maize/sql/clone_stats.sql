select count(*) from contig;

select clone.status, min(c.size), max(c.size), avg(c.size), count(*)
  from clone
       left join (select id, clone_id, end_coord-start_coord+1 as size from contig) c
       on c.clone_id = clone.id
 group by clone.status with rollup

select status, min(cnt), max(cnt), avg(cnt), sum(cnt)
  from (
select clone.status, count(*) cnt
  from clone left join contig
   on contig.clone_id = clone.id
 group by clone.id) counts
 group by status with rollup


