#!/bin/sh


#######ensembl ftp
#ensembl_archive_xx, 

######ensemblgenome ftp
#ensemblgenomes_info_xx, now ensembl_metadata 
#ensemblgenomes_stable_ids_xx_xx,
#ensembl_compara_pan_homology_xx_xx,
#ensembl_ontology_xx,
#ensembl_production_xx
#ensembl_website_xx
#ontology_mart_xx
#ncbi_taxonomy


#ftp://ftp.ensemblgenomes.org/pub/release-39/pan_ensembl/mysql/
export EG_FTP="ftp://ftp.ensemblgenomes.org/pub/"

echo "Set ensembl version e (such as 91)"
read e
echo "set e to be $e"

echo "Set ensemblgenome version eg (such as 38)"
read eg
echo "set eg to be $eg"

echo "Set gramene version g (such as 57)"
read g
echo "set g to be $g"

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
	#for db in  ensembl_compara_pan_homology_${eg}_$e ensemblgenomes_stable_ids_${eg}_$e ensembl_metadata ensembl_ontology_$e ensembl_production_$e ensembl_website_$e ontology_mart_$e 
	#	do echo "donwload $db"
	#	wget -r -nH --cut-dir=4 ${EG_FTP}/release-${eg}/pan_ensembl/mysql/$db
	#done
	wget -r -nH --cut-dir=4 ${EG_FTP}/release-${eg}/pan_ensembl/mysql/*
	~/scripts/phh-file-rename s/${eg}/${g}/ mysql/* 
	cd ../
fi


