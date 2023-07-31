#!/bin/bash

prefix=$1
gene=$2
dirr=$3
max_jobs=$4
user=$5

cd ${dirr}/${gene}_out

echo "sample	sequence	reads	identity	taxa	taxid	phylum	class	order	family	genus	bitscore	tax_num	all_species_in_best_hit" >  ${prefix}_${gene}_ttb_header
ls *_seqs.txt > samplist
num_seqs=$( cat samplist | wc -l | awk '{print $1}')
tot_per_file=$(( $num_seqs / $max_jobs ))
if [[ ${tot_per_file} -eq 0 ]]
then
  tot_per_file=1
fi
echo "there were $num_seqs samples to do taxatables for and $tot_per_file sample(s) per job."
x=1
while [[ $x -lt ${max_jobs} ]];
do
	if [[ -s samplist ]];
 	then
  		head -n ${tot_per_file} samplist > samplist_${x}
    		sed -i "1,${tot_per_file}d" samplist
      		x=$(( $x + 1 ))	
	else
 		x=$max_jobs
  	fi
done
rm samplist 	
while true;
do
	for fil in samplist_*
	do
		echo "making taxtable, doing $fil"
  		y=$( echo $ fil | awk -F"_" '{print $2}')
		if [[ ! -f ${prefix}_${gene}_taxatable.txt_${y} ]] | [[ ! -s ${prefix}_${gene}_taxatable.txt_${y} ]] ;
  		then
  			echo "outfile for $fil does not yet exist or is empty. Doing $fil."
     			sbatch ${dirr}/scripts/run_ttb.sh $fil $dirr $gene $prefix
		fi
	done
 	while true;
	do
       		sleep 2s
       	 	ck="squeue -u ${user}"
        	chck=$($ck)
       		check=$(echo "$chck" | grep "ttb" | wc -l | awk '{print $1}')
       		echo "there are $check jobs left"
		if [[ $check -eq 0 ]];then
           		echo "done with jobs, checking all ran" 
           		break
        	fi 
	done
  	ls *_taxatable.txt_* > outs_${gene}
   	num_outs=$( wc -l outs_${gene} | awk '{print $1}')
    	ls samplist_* > outs_samps
     	num_samps=$( wc -l outs_samps | awk '{print $1}')
      	if [[ $num_outs -eq $num_samps ]]; then
       		echo "all outfiles created and not empty."
	 	break
	fi
done



cat ${prefix}_${gene}_taxatable.txt_* > ttb
cat ${prefix}_${gene}_ttb_header ttb > ${prefix}_${gene}_unfiltered_taxatable.txt
rm ttb ${prefix}_${gene}_ttb_header ${prefix}_${gene}_taxatable.txt_* samplist_*

num_samp2=$( cut -f1 ${prefix}_${gene}_unfiltered_taxatable.txt | sort | uniq | wc -l | awk '{print $1}')
if [[ $num_samp2 -eq $num_seqs ]]; then
	echo "all seqfiles successfully incorporated into taxatable."
else
 	echo "there are $num_seqs seqfiles but only $num_samp2 samples in the taxatable. Some of your jobs may not have run correctly."
fi

echo "done with unfiltered analysis ${gene}. You can find taxa tables and raw ASV tables in your project directory/reports, and if you set filter to TRUE, the next step will filter your data at your read and identity cutoffs."

rm temp*
mkdir -p ${dirr}/${gene}_out/sample_seqfiles ${dirr}/results_tables
mv *_seqs.txt ${dirr}/${gene}_out/sample_seqfiles
cd ${dirr}
cp ${dirr}/${gene}_out/*_taxatable.txt ${dirr}/results_tables/
