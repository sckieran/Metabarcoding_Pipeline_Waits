#!/bin/bash

pattern=$1
r2_pattern=$2
dir=$3
max_jobs=$4
user=$4
gene=$5

cd ${dir}/${gene}

module load pear

ls *${pattern} > seqlist
num_seqs=$( cat seqlist | wc -l | awk '{print $1}')
tot_per_file=$(( $num_seqs / $max_jobs ))

x=1
while [[ $x -lt ${max_jobs} ]]
do
  if [[ -s seqlist ]];
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
  sbatch ${dir}/scripts/pear.sh $fil $pattern $r2_pattern $dir
done

while true;
do
        sleep 30s
        ck="squeue -u ${user}"
        chck=$($ck)
        check=$(echo "$chck" | grep "pr" | wc -l | awk '{print $1}')
        if [ "$check" = "0" ];then
           echo "done with pears" 
           break
        fi 
done
