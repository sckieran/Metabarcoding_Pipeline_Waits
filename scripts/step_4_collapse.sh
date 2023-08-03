#!/bin/bash


dir=$1
pattern=$2
r2_pattern=$3
max_jobs=$4
user=$5
gene=$6

cd ${dir}/${gene}

#clean up files#
mkdir -p unpaired paired collapsed seqfiles paired/unassembled  paired/discarded
mv *${pattern} ./unpaired/
mv *${r2_pattern} ./unpaired/
mv *.discarded.fastq ./paired/discarded/
mv *.unassembled*.fastq ./paired/unassembled/

#make list of files to collapse#
ls *_paired.assembled.fastq > pairedlist
num_seqs=$( wc -l pairedlist | awk '{print $1}')
tot_per_file=$( awk -v a1=$num_seqs -v a2=$max_jobs 'BEGIN { rounded = sprintf("%.0f", a1/a2); print rounded }')
if [[ ${tot_per_file} -eq 0 ]];
then
  tot_per_file=1
fi
echo "there were $num_seqs samples to pear and $tot_per_file sample(s) per job."

#cut into slurm jobs for faster processing#
x=1
while [[ $x -lt ${max_jobs} ]];
do
  if [[ -s pairedlist ]];
  then
    head -n ${tot_per_file} pairedlist > pairedlist_${x}
    sed -i "1,${tot_per_file}d" pairedlist
    x=$(( $x + 1 ))
  else
    x=$max_jobs
  fi
done
rm pairedlist
while [[ $num_seqs -ne $num_outs ]];
do
  for fil in pairedlist_*;
  do
    echo "doing $fil"
    sbatch ${dir}/scripts/run_collapser.sh $fil ${dir}/${gene}
done

while true;
do
        sleep 5s
        ck="squeue -u ${user}"
        chck=$($ck)
        check=$(echo "$chck" | grep "fx_col" | wc -l | awk '{print $1}')
        if [[ $check -eq 0 ]];then
           echo "done with collapsing ASVs" 
           break
        fi 
done

rm pairedlist_*
mv *_paired.assembled.fastq ./paired/
