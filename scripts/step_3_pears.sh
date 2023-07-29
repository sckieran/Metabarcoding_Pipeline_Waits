#!/bin/bash

pattern=$1
r2_pattern=$2
dir=$3
max_jobs=$4
user=$5
gene=$6

cd ${dir}/${gene}


ls *${pattern} > seqlist
num_seqs=$( cat seqlist | wc -l | awk '{print $1}')
tot_per_file=$(( $num_seqs / $max_jobs ))

if [[ ${tot_per_file} -eq 0 ]]
then
  tot_per_file=1
fi
echo "there were $num_seqs samples to pear and $tot_per_file sample(s) per job."

x=1
while [[ $x -lt ${max_jobs} ]]
do
  if [[ -s seqlist ]];
  then
    echo "seqlist not empty, making seqlist_${x}
    head -n ${tot_per_file} seqlist > seqlist_${x}
    sed -i "1,${tot_per_file}d" seqlist
    rem=$( wc -l seqlist | awk '{print $1}')
    echo "there are $rem samples remaining in the seqlist"
    x=$(( $x + 1 ))
  else
    x=${max_jobs}
  fi
done
rm seqlist

for fil in seqlist_*;
do
  #sbatch ${dir}/scripts/pear.sh $fil $pattern $r2_pattern ${dir}/${gene}
done

while true;
do
        sleep 30s
        ck="squeue -u ${user}"
        chck=$($ck)
        check=$(echo "$chck" | grep "pear" | wc -l | awk '{print $1}')
        if [ "$check" = "0" ];then
           echo "done with pears" 
           break
        fi 
done
