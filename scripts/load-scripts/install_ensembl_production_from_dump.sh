#!/bin/sh


DBDIR=$1
DBHOST=$2
db=ensembl_production

if [[ -z $DBDIR ]]
then

	echo "Need working directory."
        exit
fi

if [[ -z $DBHOST ]]
then

        echo "Need DBHOST to install databases on."
        exit
fi

echo "the working directory is $DBDIR, DBHOST is $DBHOST, continue? [Y/N]"

read reply

echo "reply is $reply"

if [[  $reply != [Yy] ]]
then 
	echo "Bye"
	exit
fi

echo "change directory to $DBDIR"
cd $DBDIR

mysqlparam=" -u $USER -p$DBPASSWD -h $DBHOST "
echo "mysqlparam is $mysqlparam"


/usr/local/mysql/bin/mysql $mysqlparam -e "drop database if exists $db; create database $db"  2>&1 >> ../panDB.log 
gunzip -c *.sql.gz | /usr/local/mysql/bin/mysql $mysqlparam $db 2>&1 >> ../panDB.log
gunzip *.txt.gz
ls *.txt | xargs /usr/local/mysql/bin/mysqlimport $mysqlparam -L -l $db 2>&1 >> ../panDB.log
cd ../ 

#set is_current to 1
for tb in `/usr/local/mysql/bin/mysql $mysqlparam $db -Nq -e "show tables like 'master\_%'"`
do echo $tb
	/usr/local/mysql/bin/mysql $mysqlparam $db -e "update $tb set is_current=1 where is_current=0" 
done

