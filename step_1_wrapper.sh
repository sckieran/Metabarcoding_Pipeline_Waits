#!/bin/bash
#SBATCH -C "ceph"
#SBATCH -e step_1_slurm.%j.err
#SBATCH -o step_1_slurm.%j.out
#SBATCH -J step_1
#SBATCH -t 1440

dirr=$PWD #full path to your_project directory. If you're running this code from inside this folder (recommended), you can use $PWD here.
db_dirr=reference_database #name (not path) to the out directory for your reference database. Default is "reference_database". Path is ${dirr}/reference_database.
genelist=$PWD/genelist #full path to your tab-delimited list of gene/primer set search terms.
taxlist=$PWD/taxlist #full path to your list of taxa. Space-delimited, single-column, header must be called "taxnames"
prefix="your_project" #name for your outfiles. Outfiles will generally be named "your_project_gene1" for each gene/primer set.
retmax=20 #maximum number of NCBI sequences to return. Recommended value is <100. 
rlib="~/Rpackages" #your R package library path


bash ${dirr}/scripts/step_1_get_seqs_for_database.sh -n ${prefix} -d ${dirr} -g ${genelist} -t ${taxlist} -r ${retmax} -h ${db_dirr} -l ${rlib}
