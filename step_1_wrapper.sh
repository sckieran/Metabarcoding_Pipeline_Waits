#!/bin/bash
#SBATCH -C "ceph"
#SBATCH -e step_1_slurm.%j.err
#SBATCH -o step_1_slurm.%j.out
#SBATCH -J step_1
#SBATCH -t 1440



dirr=$PWD
db_dirr=reference_database
genelist=$PWD/genelist
taxlist=$PWD/taxlist
prefix="rlha_test"
retmax=20


bash step_1_make_local_database.sh -n ${prefix} -d ${dirr} -g ${genelist} -t ${taxlist} -r ${retmax} -h ${db_dirr}
