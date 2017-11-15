#!/bin/bash
 
#$ -N my_script
#$ -l m_mem_free=3G
#$ -V
#$ -v PATH
#$ -S /bin/bash
#$ -cwd

# run commands and application
./script.dump_tree_id_by_species.pl -e panzea.registry -s $1 > tree_id.$1
