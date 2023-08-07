#!/bin/bash
module load R/4.2.3
dir=$1
rlib=$2
cutoff=$3
gene=$4

cd ${dir}/${gene}

for fil in *_clustered.fasta
do
  base=$(echo $fil | awk -F"_clustered.fasta" '{print $1}')
  grep ">" $fil | awk -F"-" '{print $2}' > temp_reads_${base}
  grep -v ">" $fil > temp_seqs_${base}
  paste  temp_seqs_${base} temp_reads_${base}  > ${base}_seqs.txt
  Rscript ${dir}/scripts/filter_rra.R $rlib ${dir}/${gene} ${base}_seqs.txt $base $cutoff
  mv ${base}_seqs.txt ./unfiltered_seqfiles/
  rm temp_reads_${base} temp_seqs_${base}
done

mv *_clustered.fasta ./collapsed/
mv *_seqs.txt ./seqfiles/
