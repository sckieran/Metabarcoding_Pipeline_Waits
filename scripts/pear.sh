#!/bin/bash

#SBATCH -J pear_j
#SBATCH -e pears.%j.err
#SBATCH -o pears.%j.out
#SBATCH -C "ceph"

infil=$1
pattern=$2
r2_pattern=$3
dir=$4

cd ${dir}
module load pear
while read fil;
do
  base=$(echo $fil | awk -F"$pattern" '{print $1}')
  pear -f ${fil} -r ${base}${r2_pattern} -o ${base}_paired -j 4
done < $infil
