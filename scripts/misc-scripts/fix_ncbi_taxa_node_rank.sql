select ntn.* from ncbi_taxa_node ntn join genome_db g using(taxon_id);
update ncbi_taxa_node ntn join genome_db g using(taxon_id) set ntn.rank='species';
