#!/bin/sh

#cat ../DBsColdenReload | while read FDB TDB; do echo "$FDB $TDB"; mysql -u weix -pwarelab -h colden -e "drop database $TDB; create database $TDB"; cd $FDB; mysql -u weix -pwarelab -h colden $TDB < $FDB.sql; ls *txt | xargs mysqlimport -u weix -pwarelab -h colden  -L -l $TDB; cd /scratch/weix/Grm44/EPl; done


DBlist=$1
Host=$2
WKdir=$3

if [ ! ${DBlist} ]
   then 
        echo "please pass database list"
        exit
fi

if [ ! ${Host} ]
   then
        echo "please pass database host name (for example: colden, cabot)"
        exit
fi

if [ ! ${WKdir} ]
   then
        echo "please pass working diretory path which should be the parenet directory of the database dump you are loading from"
        exit
fi



echo "Database list file is ${DBlist}"
echo "Database host is $Host"
echo "Your working directory which should be the database dump parent directory is $WKdir"

read -p "Continue? (y/n) " -n 1
echo
if [[ $REPLY =~ ^[Yy]$ ]]
then
        echo
else
        exit
fi

 cat $DBlist | while read FDB TDB
 do  
     cd $WKdir
     echo "$FDB $TDB"
     mysql -u weix -pwarelab -h $Host -e "drop database if EXISTS $TDB; create database $TDB"
     cd $FDB
     mysql -u weix -pwarelab -h $Host $TDB < $FDB.sql
     ls *txt | xargs mysqlimport -u weix -pwarelab -h $Host  -L -l $TDB
 done

