#!/bin/bash

#SBATCH -J mksq
#SBATCH -e mksq.%j.out
#SBATCH -o mksq.%j.err
#SBATCH -C "ceph"

inp=$1
dir=$2
gene=$3
cutoff=$4
minlen=$5
env_name=$6

cd ${dir}/${gene}
eval "$(conda shell.bash hook)"
conda activate $env_name
pref=$CONDA_PREFIX

while read fil;
do
  if [[ -s $fil ]];
  then
    base=$(echo $fil | awk -F"_clustered.fasta" '{print $1}')
    grep ">" $fil | awk -F"-" '{print $2}' > temp_reads_${base}
    grep -v ">" $fil  > temp_seqs_${base}
    paste  temp_seqs_${base} temp_reads_${base}  > ${base}_seqs.txt
    cat ${base}_seqs.txt | awk -v m=$minlen '{ if (length($1) > m) print }' > temp_${base}_seqs.txt
    mv temp_${base}_seqs.txt ${base}_seqs.txt
    Rscript ${dir}/scripts/filter_rra.R ${dir}/${gene} ${base}_seqs.txt $base $cutoff ${pref}/lib/R/library/
    mv ${base}_seqs.txt ./unfiltered_seqfiles/${base}_seqs.txt
    rm temp_reads_${base} temp_seqs_${base}
  else
    echo "$fil is empty or does not exist. Omitting $fil to avoid downstream errors. Check infiles and error logs for step 3 and 4 for $fil."
    echo -n "" > ${base}_filtered_seqs.txt
    echo -n "" > ${base}_seqs.txt
  fi
done < ${inp}

