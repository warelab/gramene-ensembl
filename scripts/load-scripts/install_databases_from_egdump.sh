#!/bin/sh


DBDIR=$1
DBHOST=$2

mysql="/usr/bin/mysql"
mysqlimport="/usr/bin/mysqlimport"

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


for k in core funcgen variation otherfeature compara ancestral mart 
	do echo $k
	ls -d *$k* > ../Epl$k
done

for k in core funcgen variation otherfeature compara ancestral mart 
#for k in mart 
        do echo Epl$k 
	for db in `cat ../Epl$k `
		do echo $db  >> ../Epl$k.log
		$mysql $mysqlparam -e "drop database if exists $db; create database $db"  2>&1 >> ../Epl$k.log 
		cd $db 2>&1 >> ../../Epl$k.log
		gunzip -c *.sql.gz | $mysql $mysqlparam $db 2>&1 >> ../../Epl$k.log
		gunzip *.txt.gz 2>&1 >> ../../Epl$k.log
		ls *.txt | xargs $mysqlimport $mysqlparam -L -l $db 2>&1 >> ../../Epl$k.log
		cd ../ 2>&1 >> ../../Epl$k.log
	done
done
