update gene set analysis_id=102 where analysis_id=171;
update transcript set analysis_id=102 where analysis_id=171;
update analysis_description set displayable=0 where analysis_id=171;
update analysis_description set web_data="{'label_key' => '[biotype]', 'colour_key' => '[biotype]', 'key' => 'gramene','default' => {'MultiTop' => 'gene_label','contigviewbottom' => 'transcript_label','MultiBottom' => 'collapsed_label','contigviewtop' => 'gene_label','cytoview' => 'gene_label','alignsliceviewbottom' => 'as_collapsed_label'}, 'name' => 'mRNA alignments', 'caption' => 'mRNA alignments', 'multi_name' => 'mRNA alignments'}" where analysis_id=107;
update analysis_description set web_data="{'label_key' => '[biotype]', 'colour_key' => '[biotype]', 'key' => 'gramene','default' => {'MultiTop' => 'gene_label','contigviewbottom' => 'transcript_label','MultiBottom' => 'collapsed_label','contigviewtop' => 'gene_label','cytoview' => 'gene_label','alignsliceviewbottom' => 'as_collapsed_label'}, 'name' => 'Rejected loci', 'caption' => 'Rejected loci', 'multi_name' => 'Rejected loci'}" where analysis_id=102;
update gene set biotype='low_confidence' where biotype='protein_coding' and analysis_id=102;
update transcript set biotype='low_confidence' where biotype='protein_coding' and analysis_id=102;
