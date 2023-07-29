#!/bin/bash

#SBATCH -J $1 
#SBATCH -e ${1}.err
#SBATCH -o ${1}.out

infil=$1
pattern=$2
r2_pattern=$3

while read p;
do
  base=$(echo $fil | awk -F"$pattern" '{print $1}')
  pear -f ${fil} -r ${base}${r2_pattern} -o ${base}_paired -j 4
done < $infil
