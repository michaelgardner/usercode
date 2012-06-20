#!/bin/sh

file=$1 # up to l
awk '{print $1}' $file > tmp.txt
exec < tmp.txt
set i=1
while read line
do
    awk -v p=$line -v p1=$(pwd) '{gsub("_runcfg_",p); gsub("_in_dir_",p1); print}' mjob_tmp.sh > $line.sh
    bsub -R "pool>10000" -q cmscaf1nh -J $line < $line.sh
    i=$(expr $i + 1);
    
done
rm tmp.txt
