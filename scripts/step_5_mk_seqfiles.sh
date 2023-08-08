#!/bin/bash
module load R/4.2.3
dir=$1
rlib=$2
cutoff=$3
gene=$4
max_jobs=$5

cd ${dir}/${gene}

ls *_clustered.fasta > collapselist
num_seqs=$( wc -l collapselist | awk '{print $1}')
tot_per_file=$( awk -v a1=$num_seqs -v a2=$max_jobs 'BEGIN { x+=(a1/a2); printf("%.0f", (x == int(x)) ? x : int(x)+1) }')
if [[ ${tot_per_file} -eq 0 ]];
then
  tot_per_file=1
fi
echo "there were $num_seqs samples to make seqfiles for and $tot_per_file sample(s) per job."

#cut into slurm jobs for faster processing#
x=1
while [[ $x -lt ${max_jobs} ]];
do
  if [[ -s collapselist ]];
  then
    head -n ${tot_per_file} collapselist > collapselist_${x}
    sed -i "1,${tot_per_file}d" collapselist
    x=$(( $x + 1 ))
  else
    x=$max_jobs
  fi
done
rm collapselist
while [[ $num_seqs -ne $num_outs ]];
do
  for fil in collapselist_*;
  do
    while read p;
    do
      base=$( echo $p | awk -F"_clustered.fasta" '{print $1}')
      if [[ ! -f ${base}_filtered_seqs.txt ]];
      then
        echo "$p" >> temp_$fil
      fi
    done < ${fil}
    if [[  -s temp_${fil} ]];
    then
    mv temp_${fil} ${fil}
      while true;
     	do
     		echo "outfile for $fil does not yet exist or is empty. Doing $fil."
     		res=$(sbatch ${dir}/scripts/run_seqs.sh $fil ${dir} ${gene} ${rlib} ${cutoff})
   		if squeue -u $user | grep -q "${res##* }"; 
   		then
   		  echo "job ${res##* } for $fil submitted successfully."
       			break
     		elif [[ -f ttb.${res##* }.err ]];
	  	then
	  		echo "job ${res##* } for $fil submitted successfully."
     			break
    		else
	 		echo "job ${res##* } did not submit. Trying again."
	 	fi
  	  done
    else
       echo "all samples for $fil completed already."
    fi
  done
  while true;
  do
        sleep 3s
        ck="squeue -u ${user}"
        chck=$($ck)
        check=$(echo "$chck" | grep "fx_col" | wc -l | awk '{print $1}')
        if [[ $check -eq 0 ]];then
           echo "done with collapsing ASVs" 
           break
        fi 
  done
  ls *_filtered_seqs.txt > outslist
  num_outs=$( wc -l outslist | awk '{print $1}')
done


mv *_clustered.fasta ./collapsed/
mv *_filtered_seqs.txt ./seqfiles/
