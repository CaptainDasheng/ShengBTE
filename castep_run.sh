#!/bin/bash

max_jobs=20                   # maximum number of simultaneous jobs
cpus=1                       # number of cores for mpirun -np
jns_ls=list.txt              # temporary list where jobnames are stored

if [ -z "$1" ]; then
    echo ---------------------------------------------------
    echo Argument not provided.
    echo Usage: $0 seedname
    exit
fi

# If a file with a similar name exists delete it
if [ -f $jns_ls ]; then
    rm -rf $jns_ls
fi

seedname=$1
for i in $seedname-3RD/*; do
	echo $i/$seedname >> $jns_ls  
done

cat $jns_ls | xargs -L 1 -P $max_jobs mpirun -np $cpus castep 

rm -rf $jns_ls
