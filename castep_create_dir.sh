#!/bin/bash

if [ -z "$1" ]; then
    echo ---------------------------------------------------
    echo Argument not provided.
    echo Usage: $0 seedname
    exit
fi

seedname=$1
rm -rf $seedname-3RD 
mkdir $seedname-3RD 
for i in 3RD.$seedname.*; do 
   s=$(echo $i|cut -d"." -f3) &&
   d=$seedname-3RD/job-$s &&
   mkdir $d &&
   mv $i $d/$seedname.cell &&
   cp $seedname.param $d/$seedname.param
done
