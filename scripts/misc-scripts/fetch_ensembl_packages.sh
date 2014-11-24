#!/bin/sh

if [ ! ${branch_e} ]
   then 
      	echo "please set branch_e"
	exit
fi

echo "branch_e is ${branch_e}"

read -p "Continue? (y/n) " -n 1
echo 
if [[ $REPLY =~ ^[Yy]$ ]] 
then
	echo 
else
	exit
fi

 for repo in \
     ensembl \
     ensembl-compara \
     ensembl-funcgen \
     ensembl-variation \
     ensembl-webcode \
     ensembl-orm \
     public-plugins \
     ensembl-tools;
 do
     if [ ! -d "$repo" ]; then
         echo "Checking out $repo (branch $branch_e)"
         git clone --branch ${branch_e} https://github.com/Ensembl/${repo}
     else
         echo Already got $repo, attempting to pull...
         cd $repo
         git pull
         git status
         cd ../
     fi
     
     echo
     echo
 done

