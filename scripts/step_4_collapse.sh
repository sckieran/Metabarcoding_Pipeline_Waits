#!/bin/bash
module load fastx

dir=$1

cd ${dir}

for fil in *_paired.assembled.fastq
do
  base=$( echo $fil | awk -F"_paired.assembled.fastq" '{print $1}')
  echo "doing $base"
  fastx_collapser -v -i $fil -o ${base}_clustered
  doub_num=$( grep -m1 -n ">*-2$" ${base}_clustered | awk -F":" '{print $1}')
  head_num=$(( $doub_num - 1 ))
  head -n $head_num ${base}_clustered > ${base}_clustered.fasta
  rm ${base}_clustered
done

mv *_paired.assembled.fastq ./paired/
