#!/bin/sh


for d in export  failure_dir  img  persistent  procedure  temporary  udcCache
do echo $d
   mkdir -p tmp/$d
done
