#!/bin/sh

root=$1
from=$2
db=$3

for p in `ls $root/patch_${2}*`
do echo $p
mysql -u ensembl_rw -p'()ryz@' -h cabot $db < $p
done

