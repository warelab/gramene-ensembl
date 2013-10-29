#!/bin/sh

cd /home/weix/data/OBathiiGenome/03102009/barthii_agps/

mysql -u weix -p -h cabot -e 'drop database oryza_barthii_core_52_30; create database oryza_barthii_core_52_30;'

mysql -u weix -p -h cabot oryza_barthii_core_52_30 < /usr/local/ensembl-live/ensembl/sql/table.sql

perl /usr/local/gramene/scripts/ensembl/load-scripts/load_assembly_from_fasta.pl --ensembl_registry=/usr/local/gramene/conf/ensembl.registry --species=Oryza_barthii --assembly_version=BAC_pool_2008 --coord_system=chromosome --rank=1 -not_seq_level Obarthii_v3.0_genbank.fna

perl /usr/local/gramene/scripts/ensembl/load-scripts/load_assembly_from_fasta.pl --ensembl_registry=/usr/local/gramene/conf/ensembl.registry --species=Oryza_barthii --assembly_version=BAC_pool_2008 --coord_system=superscaffold --rank=2 -not_seq_level OB_superscaffolds.fna

perl /usr/local/gramene/scripts/ensembl/load-scripts/load_assembly_from_fasta.pl --ensembl_registry=/usr/local/gramene/conf/ensembl.registry --species=Oryza_barthii --assembly_version=BAC_pool_2008 --coord_system=scaffold --rank=3 -not_seq_level all_scaffolds.fna

perl /usr/local/gramene/scripts/ensembl/load-scripts/load_assembly_from_fasta.pl --ensembl_registry=/usr/local/gramene/conf/ensembl.registry --species=Oryza_barthii --assembly_version=BAC_pool_2008 --coord_system=contig --rank=4  all_contigs.fna

perl /usr/local/gramene/scripts/ensembl/load-scripts/load_assembly_from_jgi.pl -registry_file /usr/local/gramene/conf/ensembl.registry -s Oryza_barthii -csa scaffold -csc contig all_scaffolds.agp


perl /usr/local/gramene/scripts/ensembl/load-scripts/load_assembly_from_jgi.pl -registry_file /usr/local/gramene/conf/ensembl.registry  -s Oryza_barthii -csa superscaffold -csc scaffold OB_superscaffolds.agp

perl /usr/local/gramene/scripts/ensembl/load-scripts/load_assembly_from_jgi.pl -registry_file /usr/local/gramene/conf/ensembl.registry -s Oryza_barthii -csa chromosome -csc superscaffold OB_shortarm.agp

perl /usr/local/gramene/scripts/ensembl/scripts/path_shorted.pl -host cabot -port 3306 -user gramene_web -pass gram3n3 -dbname oryza_barthii_core_52_30 -from_coord contig -to_coord superscaffold > path_shorted_superscaffold.sql


#load assembly table manually
#insert mapping path into meta table


#QC

#perl /usr/local/gramene/scripts/scripts/check_assembly.pl -registry_file /usr/local/gramene/scripts/conf/ensembl.registry -species Oryza_barthii -acs scaffold -ccs contig OB3S.scaffolds.fasta -chrwrong scaffold.asm.wrong -v > scaffold.asm.out


#perl /usr/local/gramene/scripts/scripts/check_assembly.pl -registry_file /usr/local/gramene/scripts/conf/ensembl.registry -species Oryza_barthii -acs superscaffold -ccs contig  -chrwrong superscaffold.asm.wrong -v OB3S.superscaffolds.fasta > superscaffold.asm.out


#Get a seq region from the fasta file
#perl /usr/local/gramene/scripts/scripts/get_fasta_region.pl -seq_name OB_3S_SUPERSCAFFOLD_4 -start_coord 2210498 -end_coord 2210598 OB3S.superscaffolds.fasta
