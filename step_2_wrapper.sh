#!/bin/bash
#SBATCH -C "ceph"
#SBATCH -e new_step_2_slurm.%j.err
#SBATCH -o new_step_2_slurm.%j.out
#SBATCH -J new_step_2
#SBATCH -t 1440



dirr=$PWD
db_dirr=reference_database
genelist=$PWD/genelist
prefix="rlha_test"
comb=


bash step_2_get_seqs_from_ncbi.sh -n ${prefix} -d ${dirr} -g ${genelist}  -h ${db_dirr} -c ${comb}

