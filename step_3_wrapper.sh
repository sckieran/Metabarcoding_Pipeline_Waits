#!/bin/bash
#SBATCH -o step_3_slurm.%j.out
#SBATCH -e step_3_slurm.%j.err
#SBATCH -C "ceph"
#SBATCH -t 1440
#SBATCH -J step_3

dir=$PWD
prefix=your_project
db_dir=$PWD/database/
pattern1="_R1.fastq" #default is _R1.fastq, leave blank for default
pattern2="_R2.fastq" #default is _R2.fastq, leave blank for default
genelist=$PWD/genelist
params_file="params_file"
localdat=boco_12S #only put something here if your local database is not named "yourproject_gene_reference", which is the default for step 1.

bash step_3_filter_and_blast.sh -d ${dir} -n ${prefix} -g ${genelist} -p ${pattern1} -q ${pattern2} -m ${params_file} -l ${localdat}

