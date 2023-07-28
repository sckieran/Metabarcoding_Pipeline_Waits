#!/bin/bash


dir=
prefix=
rlib=
genelist=
taxlist=
retmax=
db_dirr=
key=
R1_pattern=
R2_pattern=
extra_seqs=
asv_rra=
taxa_rra=
identity_cutoff=


echo "###"
echo "###"
echo "now doing step one - downloading sequences for your local reference database."
echo "###"

bash step_1_get_seqs_for_database.sh -n ${prefix} -t ${taxlist} -g ${genelist} -d ${dir} -r ${retmax} -h ${db_dirr} -l ${rlib} -k ${key}

echo "###"
echo "###"
echo "moving on to step 2 - making your blast database from the downloaded NCBI sequences"
echo "###"

bash step_2_make_database.sh -n ${prefix} -h ${db_dirr} -g ${genelist} -d ${dir}

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
  
  bash ${dir}/scripts/step_3_pears.sh ${R1_pattern} ${R2_pattern} ${dir}/${gene}
 
  echo "###"
  echo "###"
  echo "done with pears merging. Now collapsing ASVs with fastx-collapser"
  echo "###"
  
  bash ${dir}/scripts/step_4_collapse.sh ${dir}/${gene}

  echo "###"
  echo "###"
  echo "done collapsing samples into unique ASVs. Making per-sample ASV files."
  echo "###"

  bash ${dir}/scripts/step_5_mk_seqfiles.sh ${dir}/${gene}

  echo "###"
  echo "###"
  echo "done making seqfiles. No RRA filtering was done, functionality coming soon. Beginning the local BLAST process"
  echo "###"

  bash step_6_blast.sh

bash step_7_make_taxatable
