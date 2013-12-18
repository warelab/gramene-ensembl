

truncate identity_xref;
truncate object_xref;
truncate xref;
truncate protein_feature;

truncate translation_stable_id;
truncate translation; 


truncate exon_stable_id;
truncate exon_transcript;
truncate exon;

truncate transcript_supporting_feature;
truncate splicing_transcript_pair;
truncate unconventional_transcript_association;
truncate transcript_attrib;
truncate transcript_stable_id;
truncate transcript;  


truncate gene_attrib;
truncate gene_stable_id;
delete g.*, a.* from analysis a join gene g using (analysis_id);
      



-- Afterwards, you need to clean up analysis table by
-- delete from analysis where logic_name in ('GeneModel_RiceIndica_BGI', 'Ncoils', 'Pfam', 'PIRSF', 'Prints', 'scanprosite', 'Seg', 'Signalp', 'Smart', 'Superfamily', 'Tigrfam', 'XrefExonerateDNA', 'XrefExonerateProtein');
-- delete df.*, dt.*, a.* from density_feature df join density_type dt using (density_type_id) join analysis a using (analysis_id) where a.logic_name in ('knownGeneDensity', 'geneDensity');

