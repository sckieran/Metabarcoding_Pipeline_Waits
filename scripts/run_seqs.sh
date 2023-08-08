#!/bin/bash

#SBATCH -J mksq
#SBATCH -e mksq.%j.out
#SBATCH -o mksq.%j.err

inp=$1
dir=$2
gene=$3
rlib=$4
cutoff=$5
minlen=$6

cd ${dir}/${gene}
module load R/4.2.3

while read fil;
do
  base=$(echo $fil | awk -F"_clustered.fasta" '{print $1}')
  grep ">" $fil | awk -F"-" '{print $2}' > temp_reads_${base}
  grep -v ">" $fil | awk -v m=$minlen '{ if (length($1) > m) print }' > temp_seqs_${base}
  paste  temp_seqs_${base} temp_reads_${base}  > ${base}_seqs.txt
  Rscript ${dir}/scripts/filter_rra.R $rlib ${dir}/${gene} ${base}_seqs.txt $base $cutoff
  mv ${base}_seqs.txt ./unfiltered_seqfiles/${base}_seqs.txt
  rm temp_reads_${base} temp_seqs_${base}
done < ${inp}

