#!/bin/bash
#SBATCH -C "ceph"
#SBATCH -e step_2.%j.err
#SBATCH -o step_2.%j.out
#SBATCH -t 1440
#SBATCH -J step_2



dirr=$PWD #change this to the directory that contains your reference_database folder, and your fastqs in a folder with your gene name.
genelist=$PWD/genelist
prefix="your_project" #for ease of use, this should match the name you gave in Step 1.
pattern1="_R1.fastq"
pattern2="_R2.fastq"
num_graphs=24 #change this based on how many quality profiles you want to look at. Adding more adds computational time to this script.#

bash step_2_process_reads_in_dada.sh  -d ${dirr} -n ${prefix} -g ${genelist} -p ${pattern1} -q ${pattern2} -k ${num_graphs}
