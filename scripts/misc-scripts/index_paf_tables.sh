#!/bin/sh

db=$1
tb=$2

mysql -u weix -pwarelab -h cabot $db -e " \\
CREATE INDEX qmember_id_idx ON $tb (qmember_id); \\
CREATE INDEX hmember_id_idx ON $tb (hmember_id); \\
CREATE INDEX qgenome_db_id_idx ON $tb (qgenome_db_id); \\
CREATE INDEX hgenome_db_id_idx ON $tb (hgenome_db_id); \\
CREATE INDEX evalue_idx ON $tb (evalue); \\
CREATE INDEX score_idx ON $tb (score); "
