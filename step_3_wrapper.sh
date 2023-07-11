#!/bin/bash
#SBATCH -C "ceph"
#SBATCH -e step_3.%j.err
#SBATCH -o step_3.%j.out
#SBATCH -t 1440
#SBATCH -J step_3



dirr=$PWD #path to your project directory. If you're running this code from inside your directory (recommended), you can set this to $PWD.
genelist=$PWD/genelist #path to your genelist, a tab-separated single-line list of genes/primer sets that correspond to your data folders.
prefix="your_project" #for ease of use, this should match the name you gave in Step 1.
pattern1="_R1.fastq"
pattern2="_R2.fastq"
num_graphs=24 #change this based on how many quality profiles you want to look at. Adding more adds computational time to this script.#
rlib="~/Rpackages" #path to your R packages

bash ${dirr}/scripts/step_3_quality_check_reads.sh  -d ${dirr} -n ${prefix} -g ${genelist} -p ${pattern1} -q ${pattern2} -k ${num_graphs} -l ${rlib}
