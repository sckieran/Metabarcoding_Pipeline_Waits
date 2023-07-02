#!/bin/bash
#SBATCH -o step_4_slurm.%j.out
#SBATCH -e step_4_slurm.%j.err
#SBATCH -C "ceph"
#SBATCH -t 1440
#SBATCH -J step_4_remote

dir=$PWD
prefix=your_project
pattern1="_R1.fastq" #default is _R1.fastq, leave blank for default
pattern2="_R2.fastq" #default is _R2.fastq, leave blank for default
genelist=$PWD/genelist
params_file="params_file"
rlib="~/Rpackages"

bash step_4_remote.sh -d ${dir} -n ${prefix} -g ${genelist} -p ${pattern1} -q ${pattern2} -m ${params_file} -l ${rlib}

