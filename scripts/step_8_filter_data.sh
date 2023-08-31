#!/bin/bash

module load R/4.2.3

prefix=$1
gene=$2
cutoff=$3
ident=$4
dir=$5
env_name=$6

cd ${dir}/${gene}_out
eval "$(conda shell.bash hook)"
conda activate $env_name
Rscript ${dir}/scripts/filter_id_taxa.R ${prefix}_${gene}_unfiltered_taxatable.txt ${dir}/${gene}_out $prefix $gene $cutoff $ident ${CONDA_PREFIX}/lib/R/library/

cp ${prefix}_${gene}_filtered_taxatable.txt ${dir}/results_tables/
