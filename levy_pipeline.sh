#!/bin/bash

#SBATCH -J metab_pipe
#SBATCH -o metab_pipe.%j.out
#SBATCH -e metab_pipe.%j.err
#SBATCH --mem=14G

dir=$PWD
prefix=pipetest_boco
rlib="~/Rpackages"
genelist=$PWD/genelist
taxlist=$PWD/taxlist
retmax=20
db_dirr=reference_database
key=
R1_pattern="_R1.fastq"
R2_pattern="_R2.fastq"
max_jobs=490
extra_seqs=
filter=TRUE
taxa_rra=0.005
identity_cutoff=97
minlen=70
return_low=TRUE
user=sblair


echo "###"
echo "###"
echo "now doing step one - downloading sequences for your local reference database."
echo "###"

bash ${dir}/scripts/step_1_get_seqs_for_database.sh -n ${prefix} -t ${taxlist} -g ${genelist} -d ${dir} -r ${retmax} -h ${db_dirr} -l ${rlib} -k ${key}

echo "###"
echo "###"
echo "moving on to step 2 - making your blast database from the downloaded NCBI sequences"
echo "###"

bash ${dir}/scripts/step_2_make_database.sh -n ${prefix} -h ${db_dirr} -g ${genelist} -d ${dir} -e ${extra_seqs}

echo "###"
echo "###"
echo "done with step 2. Now looping over each gene/locus in your data to perform the following steps. Ensure your data is in dir/gene for each gene."
echo "###"

head -n1 $genelist | sed "s:\t:\n:g" > list_of_genes.txt
while read p;
do
  gene=$p

  echo "###"
  echo "###"
  echo "beginning pears F/R read mergers on 4 threads for gene $p"
  echo "###"
  
  bash ${dir}/scripts/step_3_pears.sh ${R1_pattern} ${R2_pattern} ${dir} ${max_jobs} ${user} ${gene}
 
  echo "###"
  echo "###"
  echo "done with pears merging. Now collapsing ASVs with fastx-collapser"
  echo "###"
  
  bash ${dir}/scripts/step_4_collapse.sh ${dir} ${pattern} ${r2_pattern} ${max_jobs} ${user} ${gene}

  echo "###"
  echo "###"
  echo "done collapsing samples into unique ASVs. Making per-sample ASV files."
  echo "###"

  bash ${dir}/scripts/step_5_mk_seqfiles.sh ${dir}/${gene}

  echo "###"
  echo "###"
  echo "done making seqfiles. No RRA filtering was done, functionality coming soon. Beginning the local BLAST process"
  echo "###"

  bash ${dir}/scripts/step_6_blast.sh -n ${prefix} -g ${gene} -d ${dir} -m ${minlen} -r ${db_dirr} -c ${identity_cutoff} -t ${return_low} -j ${max_jobs} -u ${user}

   echo "###"
   echo "###"
   echo "done with BLAST. Making raw (unfiltered) taxatable"
   echo "###"
   bash ${dir}/scripts/step_7_taxatable.sh ${prefix} ${gene} ${dir} ${max_jobs} ${user}
   
   echo "###"
   echo "###"
   if [[ $filter == "TRUE" ]]
   then
     echo "Done making raw taxatable. You chose to filter your data by relative read abundance and percent identity. Now filtering per-sample taxa at your relative read abundance cutoff of ${taxa_rra} and your identity cutoff of ${identity_cutoff}."
     bash ${dir}/scripts/step_8_filter_data.sh ${prefix} ${gene} ${taxa_rra} ${identity_cutoff} ${dir}
     echo "###"
     echo "all done with ${gene}. Moving on to next gene or exiting."
    fi
  done < list_of_genes.txt
