#!/bin/sh

ENSEMBLROOT=$1

if [[ ! -d $ENSEMBLROOT ]]
then
   echo "$ENSEMBLROOT is not a ensembl server root"
   exit
fi

cd $ENSEMBLROOT
ln -s ../biomart-perl biomart-perl
ln -s ../samtools samtools            
ln -s ../bioperl-live bioperl-live
ln -s /usr/local/apache2 apache2
mkdir -p ensembl-webcode/conf/packed
mkdir -p ensembl-webcode/htdocs/minified
mkdir tmp
mkdir logs
ln -s ensembl-webcode/conf conf
cd ensembl-webcode/
ln -s ../biomart-perl biomart-perl    
#ln -s /usr/local/biomart-perl biomart-perl
cd ../
sudo chown -R nobody:nobody tmp 
sudo chown -R nobody:nobody ensembl-webcode/cbuild
