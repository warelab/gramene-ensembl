#!/bin/sh

LOGIC_NAME=$1
DBNAME=$2
HOST=$3

if [[ ! $HOST ]]
then
	HOST="cabot"
fi

echo "logic_name for the repeat features to be deleted is $LOGIC_NAME"
echo "database is $DBNAME, host is $HOST"

mysql -u weix -pwarelab -h $HOST $DBNAME -e "delete r.*, rc.*, a.*, ad.* from repeat_feature r join repeat_consensus rc using (repeat_consensus_id) join analysis a using (analysis_id) left join analysis_description ad using (analysis_id) where a.logic_name = '$LOGIC_NAME' ";
      



