#!/bin/bash

DOWNLOAD='/home/weix/download/percona-toolkit-2.1.8/'

DB1=$1
DB2=$2
DB=$3
TDB=$4

DB1=${DB1:?plz pass a DB1 database to sync plz}
DB2=${DB2:?plz pass a DB2 database to sync plz} # Default to DB1

echo "# FROM DATABASE: '$DB1'"
echo "#  TO  DATABASE: '$DB2'"
echo

SDB=$3
TDB=$4

SDB=${SDB:?plz pass the db host for DB1 $DB1 to sync plz}
TDB=${TDB:?plz pass the db host for DB1 $DB2 to sync plz} # Default to SDB

echo "# FROM DATABASE: '$SDB'"
echo "#  TO  DATABASE: '$TDB'"
echo

## Moar command line foo anyone?

## Pick SOURCE and TARGET database instances (connection scripts)

 #SDB='mysql -u weix -pwarelab -h cabot'
 #SDB='mysql -u weix -pwarelab -h colden'

#TDB=mysql-staging-1-ensrw
#TDB=mysql-staging-2-ensrw


## RUN

## Get database details (list)
#SDBL=( $( ${SDB} details ) )
#TDBL=( $( ${TDB} details ) )

## Assign variables from the above 'details list'
S_HOST=$SDB; S_PORT=3306; S_USER='weix'; S_PASS='warelab'
T_HOST=$TDB; T_PORT=3306; T_USER='weix'; T_PASS='warelab'

## Simplify cli (Percona format)
Sx="h=$S_HOST,P=$S_PORT,u=$S_USER,p=$S_PASS"
Tx="h=$T_HOST,P=$T_PORT,u=$T_USER,p=$T_PASS"

## DEBUGGING
echo "# $Sx"
echo "# $Tx"
echo
#exit



## Run per table (a beta version does this in one step)

/usr/local/mysql/bin/mysql -h $S_HOST -u $S_USER -p${S_PASS} $DB1 -Ne 'SHOW TABLES' | while read -r table
do
    $DOWNLOAD/bin/pt-table-sync \
        --verbose \
        --print \
        ${Sx},D=${DB1},t=${table} \
        ${Tx},D=${DB2},t=${table} \
        | grep -Pv "^(INSERT|UPDATE|DELETE)" ## Just get summary information
    
    echo
done

# Alternative, limit to one table...
#   < <($SDB $DB1 -Ne 'SHOW TABLES' | grep -w analysis)

