#!/bin/bash
#SBATCH -C "ceph"
#SBATCH -e new_step_2_slurm.%j.err
#SBATCH -o new_step_2_slurm.%j.out
#SBATCH -J new_step_2
#SBATCH -t 1440



dirr=$PWD #path to your database, see step 1.
db_dirr=reference_database #name (not path) of your reference_database folder
genelist=$PWD/genelist #path to your tab-delimited gene/primer set list. 
prefix="your_project" #name, same as step 1.


bash ${dirr}/scripts/step_2_get_seqs_from_ncbi.sh -n ${prefix} -d ${dirr} -g ${genelist}  -h ${db_dirr} 

