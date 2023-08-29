#!/bin/bash

#SBATCH -J fx_col
#SBATCH -e fx_col.%j.err
#SBATCH -o fx_col.%j.out
#SBATCH -C "ceph"

module load fastx

infil=$1
dir=$2

cd $dir

while read p;
do
  base=$( echo $p | awk -F"_paired.assembled.fastq" '{print $1}')
  fastx_collapser -v -i $p -o ${base}_clustered #collapse ASVs
  doub_num=$( grep -m1 -n ">*-2$" ${base}_clustered | awk -F":" '{print $1}') #remove singletons and doubletons
  head_num=$(( $doub_num - 1 ))
  head -n $head_num ${base}_clustered > ${base}_clustered.fasta #rename files
  rm ${base}_clustered
done < $infil
