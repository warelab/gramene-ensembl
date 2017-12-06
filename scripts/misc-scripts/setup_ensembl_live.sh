#!/bin/sh

ENSEMBLROOT=$1
GRM_VERSION=$2

cd /usr/local/

if [[ ! $ENSEMBLROOT ]]
then
   echo "ENSEMBLROOT not defined, quit"
   exit
fi

if [[ ! ${GRM_VERSION} ]]
then
   echo "GRM_VERSION not defined, quit"
   exit
fi

if [[  -d $ENSEMBLROOT ]]
then  
   echo "/usr/local/$ENSEMBLROOT already exists, quit"
   exit
fi
 
if [[  -d ${GRM_VERSION} ]]
then
   echo "/usr/local/${GRM_VERSION} already exists, quit"
   exit
fi

sudo mkdir $ENSEMBLROOT
sudo git clone https://weix-cshl@github.com/warelab/gramene-ensembl ${GRM_VERSION} 

sudo chown -R weix:gramene $ENSEMBLROOT ${GRM_VERSION}

cd $ENSEMBLROOT

ln -s ../${GRM_VERSION} gramene-live

echo "Set enviroment variable branch_e (such as release/91)"
read branch_e
echo "set branch_e to be $branch_e"

echo "Set enviroment variable branch_eg (such as release/eg/38)"
read branch_eg
echo "set branch_eg to be $branch_eg"

gramene-live/scripts/misc-scripts/fetch_ensembl_packages.sh
gramene-live/scripts/misc-scripts/fetch_eg_packages.sh

ln -s ../biomart-perl biomart-perl
ln -s ../samtools samtools            
ln -s ../BioPerl-1.6.922 bioperl-live
ln -s /usr/local/apache2 apache2
ln -s ../tools_data tools_data
mkdir -p ensembl-webcode/conf/packed
mkdir -p ensembl-webcode/htdocs/minified
mkdir tmp
cd tmp/
mkdir udcCache procedure export persistent failure_dir temporary
mkdir -p temporary/tools/AssemblyConverter
mkdir -p temporary/tools/Blast
mkdir -p temporary/tools/VEP
mkdir  -p temporary/vcf_tabix

cd ../
mkdir logs
ln -s ensembl-webcode/conf conf
cd ensembl-webcode/
ln -s ../biomart-perl biomart-perl    
#ln -s /usr/local/biomart-perl biomart-perl
cd ../
sudo chown -R nobody:nobody tmp 
sudo chown -R nobody:nobody ensembl-webcode/cbuild
