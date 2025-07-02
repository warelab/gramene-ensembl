#!/bin/sh


export DBPASSWD='warelab'

DBDIR=$1
DBHOST=$2

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


	for db in `ls `
		do echo $db  >> ../panDB.log
		/usr/bin/mysql $mysqlparam -e "drop database if exists $db; create database $db"  2>&1 >> ../panDB.log 
		cd $db 2>&1 >> ../../panDB.log
		gunzip -c *.sql.gz | /usr/bin/mysql $mysqlparam $db 2>&1 >> ../../panDB.log
		gunzip -f *.txt.gz 2>&1 >> ../../panDB.log
		ls *.txt | xargs /usr/bin/mysqlimport $mysqlparam -L -l -r $db 2>&1 >> ../../panDB.log
		cd ../ 2>&1 >> ../../panDB.log
	done
