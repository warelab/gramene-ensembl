DBNAME=$1
DBUSER=$2
DBPASS=$3

if ! [[ $DBNAME && $DBUSER && $DBPASS ]]
then
	echo "Need DBNAME && DBUSER && DBPASS" && exit
fi

echo "$DBNAME && $DBUSER && $DBPASS"

mysql -h cabot -u $DBUSER -p$DBPASS -e "create database $DBNAME"
mysql -h cabot -u $DBUSER -p$DBPASS $DBNAME  < /usr/local/ensembl-live/ensembl/sql/table.sql
perl /usr/local/ensembl-live/ensembl-production/scripts/production_database/populate_production_db_tables.pl -h cabot -u $DBUSER -p $DBPASS -P 3306 -mh cabot -mu $DBUSER  -mP 3306 -mp $DBPASS -d $DBNAME -dp ./
