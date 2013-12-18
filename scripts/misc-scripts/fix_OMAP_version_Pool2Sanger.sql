update meta set meta_value = replace(meta_value, 'pool', 'Sanger') where meta_value like '%pool%';

update coord_system set version = replace(version, 'pool', 'Sanger');


