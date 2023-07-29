#!/bin/bash

pattern=$1
r2_pattern=$2
dir=$3

cd ${dir}

module load pear

ls *{pattern} > seqlist
num_seqs=$( cat seqlist | wc -l | awk '{print $1}')
tot_per_file=$(( $num_seqs / $max_jobs ))

x=1
while [[ $x -lt ${max_jobs} ]]
do
  if [[ -s diff seqlist ]]
  then
    head -n ${tot_per_file} seqlist > seqlist_${x}
    sed -i "1,${tot_per_file}d" seqlist
    x=$(( $x + 1 ))
  else
  x=$max_jobs
  fi
done
rm seqlist

for fil in seqlist_*;
do
  sbatch pear.sh $fil
done

mkdir -p unpaired paired collapsed seqfiles
mv *${pattern} ./unpaired/
mv *${r2_pattern} ./unpaired/
