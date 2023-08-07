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
tot_per_file=$( awk -v a1=$num_seqs -v a2=$max_jobs 'BEGIN { rounded = sprintf("%.0f", a1/a2); print rounded }')
echo "max_jobs is set to $max_jobs"
if [[ ${tot_per_file} -eq 0 ]];
then
  tot_per_file=1
fi
echo "there were $num_seqs samples to pear and $tot_per_file sample(s) per job."

x=1
while [[ $x -lt ${max_jobs} ]];
do
 # echo "x is $x and max_jobs is $max_jobs"
  if [[ -s seqlist ]];
  then
   # echo "seqlist not empty, making seqlist_${x}"
    head -n ${tot_per_file} seqlist > seqlist_${x}
    sed -i "1,${tot_per_file}d" seqlist
    rem=$( wc -l seqlist | awk '{print $1}')
   # echo "there are $rem samples remaining in the seqlist"
    x=$(( $x + 1 ))
  else
    echo "seqlist empty, moving on"
    x=${max_jobs}
  fi
done
rm seqlist
num_outs=1
while [[ $num_outs -ne $num_seqs ]];
do
  for fil in seqlist_*;
    while read p;
    do
      base=$( echo $p | awk -F"${pattern}" '{print $1}')
      if [[ ! -f ${base}_paired.assembled.fastq ]];
      then 
        echo "$p" >> temp_$fil
      fi
    done < ${fil}
    mv temp_${fil} ${fil}
    if [[ ! -s "$fil" ]]
    then
      sbatch ${dir}/scripts/pear.sh $fil $pattern $r2_pattern ${dir}/${gene}
    else
      echo "all files for $fil completed."
    fi
  done
  while true;
  do
       sleep 3s
       ck="squeue -u ${user}"
       chck=$($ck)
       check=$(echo "$chck" | grep "pear" | wc -l | awk '{print $1}')
       if [[ $check -eq 0 ]];then
          echo "done with pears" 
          break
       fi 
   done
  ls *_paired.assembled.fasta > outslist
  num_outs=$( wc -l outslist | awk '{print $1}')
done
