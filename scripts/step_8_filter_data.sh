#!/bin/bash

module load R/4.2.3

prefix=$1
gene=$2
cutoff=$3
ident=$4
dir=$5

cd ${dir}/${gene}_out

Rscript filter_id_taxa.R ${prefix}_${gene}_taxatable.txt ${dir}/${gene}_out $prefix $gene $cutoff $ident

cp ${prefix}_${gene}_filtered_taxatable.txt ${dir}/results_tables/