#!/bin/bash

#SBATCH -J metab_pipe
#SBATCH -o metab_pipe.%j.out
#SBATCH -e metab_pipe.%j.err
#SBATCH --mem=14G

dir=$PWD ##you can leave this as $PWD if you will run your pipeline_wrapper.sh script from a containing folder that contains directories for each gene containing your sequence fastq files, a directory containing all the scripts from github for the pipeline, and your taxalists and genelist. That is all you need in the folder. Otherwise, give the absolute path for that directory.
prefix=your_project #name of your project, all outfiles will start with this name. NO SPECIAL CHARACTERS ALLOWED, including "#$!*&()@^". NO SPACES ALLOWED
env_name=pipeline #name of the conda environment used to run the pipeline, see github tutorial
genelist=$PWD/genelist #name of your list of gene search terms, see example.
taxlist=$PWD/taxlist  #name of your taxalist prefix, files should be called taxalist_[gene] for each gene. You must have one for each gene, even if they are identical.
genus_search=TRUE #do you want to search at the genus level ("genus sp.") for genera in your taxa list that have no sequences for any of its species in your taxalist? 
retmax=20 #how many NCBI hits to return, recommend <100
db_dirr=reference_database  #desired name of reference database, default is "reference_database"
key=your_ncbi_key #you MUST get an NCBI key by registering for an NCBI account. See github walkthrough for details. 
R1_pattern="_R1.fastq"  #pattern that ends your forward reads. This should be consistent for all samples, [sample_name]_R1.fastq
R2_pattern="_R2.fastq"  #pattern that ends your reverse reads. This should be consistent for alls amples, [sample_name]_R2.fastq
max_jobs=10  #max number of jobs you want to submit to the cluster at one time. Unless you are comfortable with your job limits, I recomment this number is <=10.
extra_seqs="extra_seqs" #extra sequences that wouldn't be found in a genbank search that you want to add to your local database. name should be "extra_seqs_[gene]_sequences.fasta"
filter=FALSE #do you want to filter your data by taxa after taxatables are built? See walkthrough for help.
asv_rra=0.005 #do you want to pre-filter your data by ASV relative read abundance? Required for comparing local and remote results. If you don't want to pre-filter, set this value to 0.
taxa_rra=0.005 #RRA value to filter taxa by if you set filter=TRUE
identity_cutoff=97  #percent identity cutoff for taxa RRA filtering or returning BLAST hits. Use the literature for your gene/primer set to figure this out, 97-99 is common.
minlen=70 #set a minimum length for ASVs, base this value on the length of your amplicon after trimming. 
remote=FALSE #Do you want to run a remote-only BLAST search of the NCBI database?
remote_comp=FALSE #do you want to run BOTH a local and remote BLAST search?
return_low=TRUE #do you want BLAST to return hits below your identity% threshold? If TRUE, all hits will be returned. If FALSE, BLAST does not return hits below your cutoff and they will be called as "No Hit".
user=your_slurm_username #required. Usually the same as your cluster username.
email=your_ncbi_email_address #technically optional, but NCBI throws warnings if you don't include an email address##

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
