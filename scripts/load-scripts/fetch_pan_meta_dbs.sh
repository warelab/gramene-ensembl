#!/bin/sh


#######ensembl ftp
#ensembl_archive_xx, 
#ontology_mart_xx,
#ensembl_ontology_xx,
#ensembl_website_xx,

######ensemblgenome ftp
#ensemblgenomes_info_xx, 
#ensemblgenomes_stable_ids_xx_xx,
#ensembl_compara_pan_homology_xx_xx,

#ftp://ftp.ensembl.org/pub/release-92/mysql
ENSEMBL_FTP="ftp://ftp.ensembl.org/pub/"

#ftp://ftp.ensemblgenomes.org/pub/release-39/pan_ensembl/mysql/
EG_FTP="ftp://ftp.ensemblgenomes.org/pub/"

echo "Set ensembl version e (such as 91)"
read e
echo "set e to be $e"

echo "Set ensemblgenome version eg (such as 38)"
read eg
echo "set eg to be $eg"

echo "Set gramene version g (such as 57)"
read g
echo "set g to be $g"

#echo "Set ftp root site to download data (for example:ftp://ftp.ensemblgenomes.org/pub/)"
#read ftproot
#echo "set ftp root to be $ftproot"

if [[ -z $e ]]
then

        echo "Need ensembl version such as 92."
        exit
fi

if [[ -z $eg ]]
then

        echo "Need ensembl genome version such as 39."
        exit
fi

if [[ -z $g ]]
then

        echo "Need gramene version such as 57."
        exit
fi

#if [[ -z $ftproot ]]
#then

#        echo "Need ftp root to download database, such as ftp://ftp.ensemblgenomes.org/pub/."
#        exit
#fi

echo "the ensembl version is $e, the ensemlb genomes version is $eg, the gramene version is $g"
echo "The ftp url for ensembl is ${ENSEMBL_FTP}/release-${e}/mysql, continue? [Y/N]"

read reply

echo "reply is $reply"

if [[  $reply != [Yy] ]]
then
        echo "skip"
        
else
	mkdir ensembl_dbs
	cd ensembl_dbs	
	for db in ensembl_archive_$e ontology_mart_$e ensembl_ontology_$e ensembl_website_$e ensembl_production_$e ensembl_accounts ncbi_taxonomy
		do echo "download $db"
		wget -r -nH --cut-dir=3 ${ENSEMBL_FTP}/release-${e}/mysql/$db
		if [ "$db" == "ensembl_account" ]
		then
			mv mysql/ensembl_account mysql/ensembl_account_$e
		fi  
	done
	cd ../
fi 


#ftp://ftp.ensemblgenomes.org/pub/release-39/pan_ensembl/mysql/
echo "The ftp url for ensembl plants is ${EG_FTP}/release-${eg}/pan_ensembl/mysql, continue? [Y/N]"

read reply

echo "reply is $reply"

if [[  $reply != [Yy] ]]
then
        echo "skip"

else
	mkdir eg_dbs
        cd eg_dbs
	for db in  ensembl_compara_pan_homology_${eg}_$e ensemblgenomes_info_$eg ensemblgenomes_stable_ids_${eg}_$e ensembl_metadata 
		do echo "donwload $db"
		wget -r -nH --cut-dir=4 ${EG_FTP}/release-${eg}/pan_ensembl/mysql/$db
	done
	~/scripts/phh-file-rename s/${eg}/${g}/ mysql/* 
	cd ../
fi


