#!/bin/sh

dbname=$1
species=$2

mysql_cmd="mysql -u plensembl -pAudreyII -h bhsqldw2 "

$mysql_cmd $dbname -Nq -e "select concat('>', s.stable_id, '\n', seq.sequence, '\n') from seq_member s join gene_member g on (s.seq_member_id=g.canonical_member_id) join sequence seq using(sequence_id) join genome_db gd on(s.genome_db_id=gd.genome_db_id) where gd.name='$species' " >$species.input_pep.fa

