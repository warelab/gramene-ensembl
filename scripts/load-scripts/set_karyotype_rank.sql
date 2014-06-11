SELECT @krid:=attrib_type_id
FROM attrib_type WHERE code = 'karyotype_rank';
 
SELECT @i:=0;
 
## Assuming the name is simply the chromosome number...
INSERT INTO seq_region_attrib
SELECT seq_region_id, @krid, @i:=@i+1 FROM seq_region
INNER JOIN coord_system USING (coord_system_id)
WHERE coord_system.name = 'chromosome'
ORDER BY seq_region.name + 0;
 
## and to double check...
SELECT
  seq_region.name
FROM
  seq_region_attrib
INNER JOIN
  seq_region   USING (seq_region_id)
INNER JOIN
  coord_system USING (coord_system_id)
WHERE
  attrib_type_id = @krid
ORDER BY
  value
;
