#!/bin/bash
#SBATCH -o step_4_slurm.%j.out
#SBATCH -e step_4_slurm.%j.err
#SBATCH -C "ceph"
#SBATCH -t 1440
#SBATCH -J step_4

dirr=$PWD #path to the your_project directory. If you are running this code from inside that directory (recommended), you can put $PWD here.
prefix=your_project
db_dir=reference_database #name, not path, of database directory (path will be ${dir}/${db_dir}). Default is 'reference_database'
pattern1="_R1.fastq" #default is _R1.fastq, leave blank for default
pattern2="_R2.fastq" #default is _R2.fastq, leave blank for default
genelist=$PWD/genelist #path to genelist
params_file="params_file" #name of params_file prefix. All parameter filenames should start with this prefix and end with _gene1,_gene2...geneN for each gene/primer set in the genelist. Default is params_file
localdat= #only put something here if your local database is not named "yourproject_gene_reference", which is the default for steps 1 and 2.
rlib="~/Rpackages" #path to your Rpackages
cutoff=98
return_low=TRUE
score="pident"
local="local"


bash ${dirr}/scripts/step_4_start.sh -d ${dir} -n ${prefix} -g ${genelist} -p ${pattern1} -q ${pattern2} -m ${params_file} -l ${rlib} -c ${cutoff} -t ${return_low} -s ${score} -x ${local}  -b ${localdat}




