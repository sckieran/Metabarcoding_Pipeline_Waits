#!/bin/bash
#SBATCH -o step_5_slurm.%j.out
#SBATCH -e step_5_slurm.%j.err
#SBATCH -C "ceph"
#SBATCH -t 2880
#SBATCH -J step_5

dir=$PWD #path to the your_project directory. If you are running this code from inside that directory (recommended), you can put $PWD here.
prefix=your_project
genelist=$PWD/genelist #path to genelist
cutoff=97 #your identity % cutoff. Non-inclusive. 
rlib="~/Rpackages" #path to your r packages library

bash step_4_filter_and_blast_local.sh -d ${dir} -n ${prefix} -g ${genelist} -c ${cutoff} -l ${rlib}

