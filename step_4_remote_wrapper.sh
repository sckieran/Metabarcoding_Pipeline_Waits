#!/bin/bash
#SBATCH -o step_3_slurm.%j.out
#SBATCH -e step_3_slurm.%j.err
#SBATCH -C "ceph"
#SBATCH -t 1440

dir=$PWD
prefix=boco
pattern1="_R1.fastq" #default is _R1.fastq, leave blank for default
pattern2="_R2.fastq" #default is _R2.fastq, leave blank for default
genelist=$PWD/genelist
params_file="params_file"

bash step_4_remote.sh -d ${dir} -n ${prefix} -g ${genelist} -p ${pattern1} -q ${pattern2} -m ${params_file}

