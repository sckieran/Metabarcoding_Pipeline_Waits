#!/bin/bash

#SBATCH -J metab_pipe
#SBATCH -o metab_pipe.%j.out
#SBATCH -e metab_pipe.%j.err
#SBATCH --mem=14G

dir=$PWD
prefix=your_project
env_name=pipeline_env
genelist=$PWD/genelist
taxlist=$PWD/taxlist
genus_search=TRUE #do you want to search at the genus_sp. level for genera in your taxa list that have no species for a gene?
retmax=20 #how many NCBI hits to return
db_dirr=reference_database
key=your_ncbi_key 
R1_pattern="_R1.fastq"
R2_pattern="_R2.fastq"
max_jobs=10
extra_seqs="extra_seqs"
filter=TRUE
asv_rra=0.005
taxa_rra=0.005
identity_cutoff=97
minlen=70
remote=FALSE
remote_comp=FALSE
return_low=TRUE
user=your_slurm_username
email=your_ncbi_email_address #technically optional, but NCBI throws warnings if you don't include

echo "###"
echo "###"
echo "now doing step one - downloading sequences for your local reference database."
echo "###"

bash ${dir}/scripts/step_1_get_seqs_for_database.sh -n ${prefix} -t ${taxlist} -g ${genelist} -d ${dir} -r ${retmax} -h ${db_dirr} -s ${genus_search} -p ${env_name} -k ${key} -e ${email} 
exit_status=$?
if [ "${exit_status}" -ne 0 ];
then
    echo "Step 1 failed. exit ${exit_status}"
    exit 1
fi

echo "###"
echo "###"
echo "moving on to step 2 - making your blast database from the downloaded NCBI sequences"
echo "###"

bash ${dir}/scripts/step_2_make_database.sh -n ${prefix} -h ${db_dirr} -g ${genelist} -d ${dir} -f ${env_name} -e ${extra_seqs}
exit_status=$?
if [ "${exit_status}" -ne 0 ];
then
    echo "Step 2 failed. exit ${exit_status} and restart at step 2."
    exit 1
fi

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
  exit_status=$?
  if [ "${exit_status}" -ne 0 ];
  then
    echo "Step 3 failed. exit ${exit_status} and restart at step 3."
    exit 1
  fi
 
  echo "###"
  echo "###"
  echo "done with pears merging. Now collapsing ASVs with fastx-collapser"
  echo "###"
  
  bash ${dir}/scripts/step_4_collapse.sh ${dir} ${R1_pattern} ${R2_pattern} ${max_jobs} ${user} ${gene}
  exit_status=$?
  if [ "${exit_status}" -ne 0 ];
  then
    echo "Step 4 failed. exit ${exit_status} and restart at step 4."
    exit 1
  fi
  echo "###"
  echo "###"
  echo "done collapsing samples into unique ASVs. Making per-sample ASV files."
  echo "###"

  bash ${dir}/scripts/step_5_mk_seqfiles.sh ${dir} ${asv_rra} ${gene} ${max_jobs} ${user} ${minlen} ${env_name}
  exit_status=$?
  if [ "${exit_status}" -ne 0 ];
  then
    echo "Step 5 failed. exit ${exit_status} and restart at step 5."
    exit 1
  fi
  echo "###"
  echo "###"
  echo "done making seqfiles. RRA filtering ASV/sample done at your filter level of ${taxa_rra}. Beginning the BLAST process"
  echo "###"
  if [[ $remote == "TRUE" ]]
  then
    echo "building input FASTA and performing remote BLAST."
    bash ${dir}/scripts/step_6_blast_remote.sh -n ${prefix} -g ${gene} -d ${dir} -m ${minlen} -c ${identity_cutoff} -t ${return_low} -j ${max_jobs} -u ${user}
    exit_status=$?
    if [ "${exit_status}" -ne 0 ];
    then
      echo "Step 6 failed. exit ${exit_status} and restart at step 6."
      exit 1
    fi
    echo "###"
    echo "###"
    echo "done with BLAST. Making raw (unfiltered) taxatable"
    echo "###"
    bash ${dir}/scripts/step_7_taxatable.sh ${prefix} ${gene} ${dir} ${max_jobs} ${user}
    if [ "${exit_status}" -ne 0 ];
    then
      echo "Step 7 failed. exit ${exit_status} and restart at step 7."
      exit 1
    fi
elif [[ $remote_comp == "TRUE" ]]
    then
    echo "performing both local and remote BLAST and comparing the two taxatables."
    bash ${dir}/scripts/step_6_blast_comparison.sh -n ${prefix} -g ${gene} -d ${dir} -m ${minlen} -r ${db_dirr} -c ${identity_cutoff} -t ${return_low} -j ${max_jobs} -u ${user} -e ${env_name}
    exit_status=$?
    if [ "${exit_status}" -ne 0 ];
    then
        echo "Step 6 failed. exit ${exit_status} and restart at step 6."
        exit 1
    fi
    echo "###"
    echo "###"
    echo "done with BLAST. Making raw (unfiltered) taxatable"
    echo "###"
    bash ${dir}/scripts/step_7_taxatable_comp.sh ${prefix} ${gene} ${dir} ${max_jobs} ${user}
    if [ "${exit_status}" -ne 0 ];
    then
      echo "Step 7 failed. exit ${exit_status} and restart at step 7."
      exit 1
    fi
  else
    echo "building input FASTA and performing local BLAST."
    bash ${dir}/scripts/step_6_blast.sh -n ${prefix} -g ${gene} -d ${dir} -m ${minlen} -r ${db_dirr} -c ${identity_cutoff} -t ${return_low} -j ${max_jobs} -u ${user}
    exit_status=$?
    if [ "${exit_status}" -ne 0 ];
    then
      echo "Step 6 failed. exit ${exit_status} and restart at step 6."
      exit 1
    fi
   echo "###"
   echo "###"
   echo "done with BLAST. Making raw (unfiltered) taxatable"
   echo "###"
   bash ${dir}/scripts/step_7_taxatable.sh ${prefix} ${gene} ${dir} ${max_jobs} ${user}
    if [ "${exit_status}" -ne 0 ];
    then
      echo "Step 7 failed. exit ${exit_status} and restart at step 7."
      exit 1
    fi  
  fi

   
   if [[ $filter == "TRUE" ]]
   then
     echo "Done making raw taxatable. You chose to filter your data by relative read abundance and percent identity. Now filtering per-sample taxa at your relative read abundance cutoff of ${taxa_rra} and your identity cutoff of ${identity_cutoff}."
     if [[ $remote_comp == "TRUE" ]]
     then
         echo "This option cannot be used with remote_compare=TRUE. Skipping this step."
         exit 1
     fi
     bash ${dir}/scripts/step_8_filter_data.sh ${prefix} ${gene} ${taxa_rra} ${identity_cutoff} ${dir} ${env_name}
     exit_status=$?
     if [ "${exit_status}" -ne 0 ];
     then
        echo "Step 5 failed. exit ${exit_status} and restart at step 5."
        exit 1
     fi
     echo "###"
     echo "all done with ${gene}. Moving on to next gene or exiting."
    fi
  done < list_of_genes.txt
