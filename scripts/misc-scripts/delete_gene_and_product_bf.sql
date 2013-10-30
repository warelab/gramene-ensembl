
truncate dependent_xref;
truncate ontology_xref;
truncate identity_xref;
truncate object_xref;
truncate xref;
truncate protein_feature;

truncate translation; 


truncate exon_transcript;
truncate exon;

truncate transcript_supporting_feature;
truncate splicing_transcript_pair;
truncate unconventional_transcript_association;
truncate transcript_attrib;
truncate transcript;  


truncate gene_attrib;
delete g.*, a.*, ad.* from analysis a left join analysis_description ad using (analysis_id) join gene g using (analysis_id);
      



-- Afterwards, you need to clean up analysis and analysis_description tables by
-- delete a.*, ad.* from analysis a join analysis_description ad using (analysis_id) where logic_name in ('GeneModel_RiceIndica_BGI', 'Ncoils', 'Pfam', 'PIRSF', 'Prints', 'scanprosite', 'Seg', 'Signalp', 'Smart', 'Superfamily', 'Tigrfam', 'XrefExonerateDNA', 'XrefExonerateProtein');
-- delete df.*, dt.*, a.* from density_feature df join density_type dt using (density_type_id) join analysis a using (analysis_id) where a.logic_name in ('knownGeneDensity', 'geneDensity');

