#!/bin/sh

if [ ! ${branch_e} ]
   then 
      	echo "please set branch_e, such as release/79"
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
     ensembl-production \
     ensembl-compara \
     ensembl-funcgen \
     ensembl-variation \
     ensembl-webcode \
     ensembl-orm \
     public-plugins \
     ensembl-tools \
     ensembl-hive \
     ensembl-production \
     ensembl-io \
<<<<<<< HEAD
     ensembl-rest
     ensembl-variation/VEP_plugins;
=======
     VEP_plugins;
>>>>>>> b7c6b6135e50605224278b077411719a1d42537e
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

