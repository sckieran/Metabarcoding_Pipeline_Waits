#!/bin/bash

dir=$1

cd ${dir}

for fil in *_clustered.fasta
do
  base=$(echo $fil | awk -F"_clustered.fasta" '{print $1}')
  grep ">" $fil | awk -F"-" '{print $2}' > temp_reads_${base}
  grep -v ">" $fil > temp_seqs_${base}
  paste  temp_seqs_${base} temp_reads_${base}  > ${base}_seqs.txt
  #Rscript filter_rra.R ${base}_seqs.txt
  rm temp_reads_${base} temp_seqs_${base}
done

mv *_clustered.fasta ./collapsed/
