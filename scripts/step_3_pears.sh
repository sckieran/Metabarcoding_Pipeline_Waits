#!/bin/bash

pattern=$1
r2_pattern=$2
dir=$3

cd ${dir}

module load pear

for fil in *_${pattern};
do
  base=$(echo $fil | awk -F"$pattern" '{print $1}')
  pear -f ${fil} -r ${base}_${r2_pattern} -o ${base}_paired -j 4
done

mkdir -p unpaired paired collapsed seqfiles
mv *_${pattern} ./unpaired/
mv *_${r2_pattern} ./unpaired/
