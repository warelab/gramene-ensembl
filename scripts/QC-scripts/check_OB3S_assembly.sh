#!/bin/sh

cd /home/weix/data/OBathiiGenome/03102009/barthii_agps

perl /usr/local/gramene/scripts/ensembl/QC-scripts/check_assembly.pl -registry_file /usr/local/gramene/conf/ensembl.registry -species Oryza_barthii -acs scaffold -ccs contig all_scaffolds.fna -chrwrong scaffold.asm.wrong -v > scaffold.asm.out

perl /usr/local/gramene/scripts/ensembl/QC-scripts/check_assembly.pl -registry_file /usr/local/gramene/conf/ensembl.registry -species Oryza_barthii -acs superscaffold -ccs contig OB_superscaffolds.fna -chrwrong superscaffold.asm.wrong -v > superscaffold.asm.out

perl /usr/local/gramene/scripts/ensembl/QC-scripts/check_assembly.pl -registry_file /usr/local/gramene/conf/ensembl.registry -species Oryza_barthii -acs chromosome -ccs contig Obarthii_v3.0_genbank.fna -chrwrong chr.asm.wrong -v > chr.asm.out
                                                                          

