#!/bin/sh

if [ ! ${branch_eg} ]
   then 
      	echo "please set branch_eg, for example release/eg/26"
	exit
fi

echo "branch_eg is ${branch_eg}"
read -p "Continue? (y/n) " -n 1
echo
if [[ $REPLY =~ ^[Yy]$ ]]
then   
        echo
else   
        exit
fi


for repo in \
     eg-web-search \
     eg-web-common \
     ensemblgenomes-api \
     eg-web-plants ;
 do
     if [ ! -d "$repo" ]; then
         echo "Checking out $repo (branch $branch_eg)"
         git clone --branch ${branch_eg} https://github.com/EnsemblGenomes/${repo}
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
