
DELIMITER $$
DROP PROCEDURE IF EXISTS `truncate_if_exist`$$

CREATE  PROCEDURE `truncate_if_exist`(IN tbl_name VARCHAR(150) )
  BEGIN
    IF EXISTS( SELECT 1 FROM information_schema.TABLES WHERE table_name = tbl_name AND table_schema = DATABASE()) THEN
    SET @query = CONCAT('TRUNCATE ', tbl_name);
    PREPARE stmt FROM @query;
    EXECUTE stmt;
    DEALLOCATE PREPARE stmt;
    END IF;
  END $$

DELIMITER ;

CALL truncate_if_exist('dependent_xref');
CALL truncate_if_exist('ontology_xref');
CALL truncate_if_exist('identity_xref');
CALL truncate_if_exist('object_xref');
CALL truncate_if_exist('xref');
CALL truncate_if_exist('protein_feature');
CALL truncate_if_exist('translation');
CALL truncate_if_exist('exon_transcript');
CALL truncate_if_exist('exon');
CALL truncate_if_exist('transcript_supporting_feature');
CALL truncate_if_exist('splicing_transcript_pair');
CALL truncate_if_exist('transcript_attrib');
CALL truncate_if_exist('transcript');
CALL truncate_if_exist('gene_attrib');

-- truncate dependent_xref;
-- truncate ontology_xref;
-- truncate identity_xref;
-- truncate object_xref;
-- truncate xref;
-- truncate protein_feature;

-- truncate translation; 


-- truncate exon_transcript;
-- truncate exon;

-- truncate transcript_supporting_feature;
-- truncate splicing_transcript_pair;
-- truncate transcript_attrib;
-- truncate transcript;  


-- truncate gene_attrib;

delete g.*, a.*, ad.* from analysis a left join analysis_description ad using (analysis_id) join gene g using (analysis_id);
      



-- Afterwards, you need to clean up analysis and analysis_description tables by
-- delete a.*, ad.* from analysis a join analysis_description ad using (analysis_id) where logic_name in ('GeneModel_RiceIndica_BGI', 'Ncoils', 'Pfam', 'PIRSF', 'Prints', 'scanprosite', 'Seg', 'Signalp', 'Smart', 'Superfamily', 'Tigrfam', 'XrefExonerateDNA', 'XrefExonerateProtein');
-- delete df.*, dt.*, a.* from density_feature df join density_type dt using (density_type_id) join analysis a using (analysis_id) where a.logic_name in ('knownGeneDensity', 'geneDensity');

