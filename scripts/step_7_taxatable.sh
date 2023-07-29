#!/bin/bash

prefix=$1
gene=$2
dirr=$3
max_jobs=$4
user=$5

cd ${dirr}/${gene}_out

echo "sample	sequence	reads	identity	taxa	taxid	phylum	class	order	family	genus	bitscore	tax_num	all_species_in_best_hit	seqnum" >  ${prefix}_${gene}_ttb_header
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

for fil in samplist_*
do
	echo "making taxtable, doing samplist_${x}"
	sbatch ${dirr}/scripts/run_ttb.sh samplist_${x} $dirr $gene $prefix
done

while true;
do
        sleep 10s
        ck="squeue -u ${user}"
        chck=$($ck)
        check=$(echo "$chck" | grep "ttb" | wc -l | awk '{print $1}')
        if [[ $check -eq 0 ]];then
           echo "done with taxatable" 
           break
        fi 
done

cat ${prefix}_${gene}_taxatable.txt_* > ttb
cat ${prefix}_${gene}_ttb_header ttb > ${prefix}_${gene}_unfiltered_taxatable.txt
rm ttb ${prefix}_${gene}_ttb_header ${prefix}_${gene}_taxatable.txt_* samplist_*

echo "done with unfiltered analysis ${gene}. You can find taxa tables and raw ASV tables in your project directory/reports, and if you set filter to TRUE, the next step will filter your data at your read and identity cutoffs."
rm temp*
mkdir ${dirr}/${gene}_out/sample_seqfiles
mv *_seqs.txt ${dirr}/${gene}_out/sample_seqfiles
cd ${dirr}
cp ${dirr}/${gene}_out/*_taxatable.txt ${dirr}/results_tables/
