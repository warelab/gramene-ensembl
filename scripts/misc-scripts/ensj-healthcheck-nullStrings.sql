UPDATE analysis SET db_version = NULL WHERE db_version = '' OR db_version = 'NULL';
UPDATE analysis SET db_file = NULL WHERE db_file = '' OR db_file = 'NULL';
UPDATE analysis SET program_version = NULL WHERE program_version = '' OR program_version = 'NULL';
UPDATE analysis SET program_file = NULL WHERE program_file = '' OR program_file = 'NULL';
UPDATE analysis SET parameters = NULL WHERE parameters = '' OR parameters = 'NULL';
UPDATE analysis SET module = NULL WHERE module = '' OR module = 'NULL';
UPDATE analysis SET module_version = NULL WHERE module_version = '' OR module_version = 'NULL';
UPDATE analysis SET gff_source = NULL WHERE gff_source = '' OR gff_source = 'NULL';
UPDATE analysis SET gff_feature = NULL WHERE gff_feature = '' OR gff_feature = 'NULL';
UPDATE analysis_description SET web_data = NULL WHERE web_data = '' OR web_data = 'NULL';
UPDATE object_xref SET linkage_annotation = NULL WHERE linkage_annotation = '' OR linkage_annotation = 'NULL';
UPDATE protein_feature SET external_data = NULL WHERE external_data = '' OR external_data = 'NULL';
UPDATE protein_feature SET hit_description = NULL WHERE hit_description = '' OR hit_description = 'NULL';
