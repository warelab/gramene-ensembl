# Determines all top-level seq_regions in the database and flags these
# as toplevel in the seq_region_attrib table

INSERT IGNORE INTO attrib_type ( code, name, description ) 
VALUES ( 'toplevel','Top Level','Top Level Non-Redundant Sequence Region' );

# Remove old
DELETE sra FROM seq_region_attrib sra
JOIN   attrib_type at using(attrib_type_id) 
WHERE  at.code='toplevel';

# Add new
REPLACE INTO seq_region_attrib
SELECT distinct( sr.seq_region_id ), at.attrib_type_id, 1
FROM   seq_region sr LEFT JOIN assembly a 
       ON sr.seq_region_id=a.cmp_seq_region_id,
       attrib_type at
WHERE  cmp_seq_region_id is NULL
AND    at.code='toplevel';

# Update assembly.num_toplevel_seqs value in meta table
DELETE FROM meta 
WHERE meta_key='assembly.num_toplevel_seqs'; 
INSERT INTO meta(meta_key,meta_value) 
SELECT 'assembly.num_toplevel_seqs', count(*) 
FROM   seq_region join seq_region_attrib using(seq_region_id) 
       join attrib_type using(attrib_type_id) 
WHERE  code='toplevel';

