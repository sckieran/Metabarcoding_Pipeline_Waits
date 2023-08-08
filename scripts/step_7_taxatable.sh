#!/bin/bash

prefix=$1
gene=$2
dirr=$3
max_jobs=$4
user=$5

cd ${dirr}/${gene}_out

echo "sample	sequence	reads	identity	taxa	taxid	phylum	class	order	family	genus	bitscore	tax_num	all_species_in_best_hit" >  ${prefix}_${gene}_ttb_header
ls *_seqs.txt > samplist
num_seqs=$( wc -l samplist | awk '{print $1}')
tot_per_file=$( awk -v a1=$num_seqs -v a2=$max_jobs 'BEGIN { x+=(a1/a2); printf("%.0f", (x == int(x)) ? x : int(x)+1) }')
if [[ ${tot_per_file} -eq 0 ]]
then
  tot_per_file=1
fi
echo "there were $num_seqs samples to do taxatables for and $tot_per_file sample(s) per job."
x=1
while [[ $x -le ${max_jobs} ]];
do
	if [[ -s samplist ]];
 	then
  		head -n ${tot_per_file} samplist > samplist_${x}
    		sed -i "1,${tot_per_file}d" samplist
      		x=$(( $x + 1 ))	
	else
 		x=$(( $max_jobs + 1 ))
  	fi
done
rm samplist 	
num_outs=1
cat *_seqs.txt | cut -f1 | awk -v m=$minlen '{ if (length($1) > m) print }' > all_outs
sed -i '/^#/d' all_outs
num_samps=$( wc -l all_outs | awk '{print $1}')
while [[ $num_outs -ne $num_samps ]];
do
 	echo "$num_outs and $num_samps are not equal, running job submission"
  	for fil in samplist_*
	do
  		y=$( echo $fil | awk -F"_" '{print $2}')
		if [[ ! -s ${prefix}_${gene}_taxatable.txt_${y} ]];
  		then
     			while true;
     			do
     				echo "outfile for $fil does not yet exist or is empty. Doing $fil."
     				res=$(sbatch ${dirr}/scripts/run_ttb.sh $fil $dirr $gene $prefix)
   				if squeue -u $user | grep -q "${res##* }"; 
   				then
   					echo "job ${res##* } for $fil submitted successfully."
       					break
       				elif [[ -f ttb.${res##* }.err ]];
	  			then
	  				echo "job ${res##* } for $fil submitted successfully."
      					break
      				else
	  				echo "job ${res##* } did not submit. Trying again."
				fi
   			done
  		else
    			touch temp_samplist_${y}
       			while read p;
      			do
				base=$( echo "$p" | awk -F"_filtered_seqs.txt" '{print $1}' )
   				num_ttb=$( grep "${base}" ${prefix}_${gene}_taxatable.txt_${y} | wc -l | awk '{print $1}' )
      				num_s=$( wc -l $p | awk '{print $1}' )
	  			echo "$base has $num_s sequences but only $num_ttb were integrated into the outfile. Resubmitting."
	 			if [[ $num_ttb -ne $num_s ]];
    				then
       					echo "$base has $num_s sequences but only $num_ttb were integrated into the outfile. Resubmitting."
	    				echo "$p" >> temp_samplist_${y}
	   			fi
      			done < samplist_${y}
      			if [[ -s temp_samplist_${y} ]];
			then
  				lensl=$( wc -l temp_samplist_${y} | awk '{print $1}' )
      				echo "There is an outfile for samplist_${y}, but it is incomplete. There are ${lensl} samples to integrate."
      				mv temp_samplist_${y} samplist_${y}
      				res=$(sbatch ${dirr}/scripts/run_ttb.sh samplist_${y} $dirr $gene $prefix)
   				if squeue -u $user | grep -q "${res##* }"; 
   				then
   					echo "job ${res##* } for $fil submitted successfully."
       					break
       				elif [[ -f ttb.${res##* }.err ]];
	  			then
	  				echo "job ${res##* } for $fil submitted successfully."
      					break
      				else
	  				echo "job ${res##* } did not submit. Trying again."
				fi
    			else
       				echo "all sequences for samplist_${y} completed."
       			fi
       		fi
	done
 	while true;
	do
       		sleep 2s
       		check=$( squeue -u $user | grep "ttb" | wc -l )
       		echo "Waiting for jobs to finish. There are $check jobs left"
		if [[ $check -eq 0 ]];then
           		echo "done with jobs, checking all ran" 
           		break
        	fi 
	done
  	cat ${prefix}_${gene}_taxatable.txt_* > outs_${gene}
   	num_outs=$( wc -l outs_${gene} | awk '{print $1}')
    	if [[ $num_outs -ne $num_samps ]]; then
     		echo "number of outfiles is $num_outs, and is not equal to the number of input sequences, $num_samps. Checking outfiles and re-running jobs."
       	else
		echo "number of outfiles, $num_outs equals the number of input samples, $num_samps. Continuing."
	fi
done



cat ${prefix}_${gene}_taxatable.txt_* > ttb
cat ${prefix}_${gene}_ttb_header ttb > ${prefix}_${gene}_unfiltered_taxatable.txt
#rm ttb ${prefix}_${gene}_ttb_header ${prefix}_${gene}_taxatable.txt_* samplist_* ttb.*

num_samp2=$( cut -f1 ${prefix}_${gene}_unfiltered_taxatable.txt | tail -n+2 | sort | uniq | wc -l | awk '{print $1}')

if [[ $num_samp2 -eq $num_seqs ]]; then
	echo "all seqfiles successfully incorporated into taxatable."
else
 	echo "there are $num_seqs seqfiles but only $num_samp2 samples in the taxatable. Some of your jobs may not have run correctly."
fi

echo "done with unfiltered analysis ${gene}. You can find taxa tables and raw ASV tables in your project directory/reports, and if you set filter to TRUE, the next step will filter your data at your read and identity cutoffs."

rm *_taxatable.txt_*
rm temp* rm ttb*.err ttb*.out samplist_*
mkdir -p ${dirr}/${gene}_out/sample_seqfiles ${dirr}/results_tables
mv *_seqs.txt ${dirr}/${gene}_out/sample_seqfiles
cd ${dirr}
cp ${dirr}/${gene}_out/*_taxatable.txt ${dirr}/results_tables/
