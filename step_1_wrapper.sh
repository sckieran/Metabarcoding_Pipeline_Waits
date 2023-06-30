#!/bin/bash
#SBATCH -C "ceph"
#SBATCH -e step_1.%j.err
#SBATCH -o step_1.%j.out
#SBATCH -J step_1
#SBATCH -t 1440



dirr=$PWD/database
genelist=$PWD/genelist
taxlist=$PWD/taxlist
prefix="boco"

bash step_1_make_local_database.sh -n ${prefix} -d ${dirr} -g ${genelist} -t ${taxlist}
