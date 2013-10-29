
delete from karyotype 
where seq_region_id in (select asm_seq_region_id from assembly);

-- Don't have true band pattern, so using clones instead
insert into karyotype 
       (seq_region_id, seq_region_start, seq_region_end, band, stain )
select sra.seq_region_id, a.asm_start, a.asm_end, src.name, 'gneg' 
from   seq_region sra, assembly a, seq_region src
where  a.asm_seq_region_id = sra.seq_region_id
and    a.cmp_seq_region_id = src.seq_region_id
order by sra.seq_region_id, a.asm_start;

-- Alter band colour for alternating clones
update karyotype
set    stain='gpos25'
where  mod( karyotype_id, 2 ) > 0.5;

-- Centromeres will be added by gramene_ensembl/tigr-fetch/centromeres.pl

-- Centromeres (approx)
-- Chr1    AP007228        B1061G08        16.8
-- Chr2    AP007145        B1120G10d       13.7
-- Chr3    AC137925        OSJNBb0047D08   19.5
-- Chr4    BX890594        OSJNBb0062N22   9.7
-- Chr5    AC137984        P0697B04        12.4
-- Chr6    AP005763        OSJNBa0015G09   15.4
-- Chr7    AP007259        B1155H11d       12.2
-- Chr8    AY360388        OSJNBa0038J12   12.9
-- Chr9    AP007147        B1106C08d       2.7
-- Chr10   AC022352        OSJNBa0034E23   7.7
-- Chr11   AC137922        OSJNBa0046A04   12.0
-- Chr12   BX000556        OSJNBa0088J04   12.0
-- 'acen'

-- Centromeres (real)
-- BA000044 (chromosome 8, centromeric region) 
-- AP005305  .3:1..177529,
-- AP005843  .2:123169..137779,
-- AP005819  .2:14138..188058,
-- AP004458  .2:73216..160541,
-- AP004041  .2:5922..165126,
-- AP005540  .2:66111..120667,
-- AP005162  .2:1624..152702,
-- AP005166  .2:110729..184235,
-- AP005693  .2:40759..155085,
-- AP006480  .2:91548..190806,
-- AP006481  .1:34347..139217,
-- AP006482  .1:35049..158084,
-- AP005832  .2:33203..176635,
-- AP004228  .2:61894..118763,
-- AP005500  .2:20351..136473,
-- AP005407  .2:55606..152373,
-- AP005498  .2:8997..163932,
-- AP004561  .2:57583..148590


select sr.name, band from karyotype k, seq_region sr 
where  sr.seq_region_id=k.seq_region_id 
and    band in ( 'AP007228','AP007145','AC137925',
                 'BX890594','AC137984','AP005763',
                 'AP007259','AY360388','AP007147',
                 'AC022352','AC137922','BX000556' );

select sr.name, band from karyotype k, seq_region sr
where  sr.seq_region_id=k.seq_region_id
and    band in ( 'AP005305','AP005843','AP005819',
                 'AP004458','AP004041','AP005540',
                 'AP005162','AP005166','AP005693',
                 'AP006480','AP006481','AP006482', 
                 'AP005832','AP004228','AP005500',
                 'AP005407','AP005498','AP004561' );

