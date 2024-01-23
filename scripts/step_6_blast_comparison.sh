#!/bin/bash

while getopts ":n:g:d:m:r:b:c:t:j:u:e:" opt; do
  case $opt in
    n) prefix="$OPTARG"
    ;;
    g) gene="$OPTARG"
    ;;
    d) dirr="$OPTARG"
    ;;
    m) minlen="$OPTARG"
    ;;
    r) db_dirr="$OPTARG"
    ;;
    b) localdat="$OPTARG"
    ;;
    c) cutoff="$OPTARG"
    ;;
    t) return_low="$OPTARG"
    ;;
    j) max_jobs="$OPTARG"
    ;;
    u) user="$OPTARG"
    ;;
    e) env_name="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    exit 1
    ;;
    esac
done


module load ncbi-blast/2.10.1

max_jobs=$(( $max_jobs / 2 )) #because we submit two jobs per set of samples, one local and one remote, we reduce the max_jobs by half so we don't go over the max. This is imperfect because of floating point, will fix later.

cd ${dirr}
#assign defaults if parameters not provided##
if [[ -z ${localdat} ]]
then
	localdat=${prefix}_${gene}_reference
fi

#make out directory#
mkdir -p ${gene}_out 

##copy filtered seqfiles, the output of step 5, to the out directory##
cd ${gene}_out
cp ${dirr}/${gene}/seqfiles/* .
 	
echo "copied files. Beginning blast and taxtable."
##get unique list of all ASVs by concatenating all seqfiles and pulling uniques from that list. Remove empty lines, this is just a formatting thing to ensure no errors##
cat *_seqs.txt | cut -f1 | sort | uniq > temp_seqs
sed -i '/^$/d' temp_seqs


##make query fasta from seqlist#
x=1
n=$(wc -l temp_seqs | awk '{print $1}')
## this loop just makes fasta headers, which start >seq_[01-100000]. It adds 0s to the front as needed and goes up to 999,999. This keeps sequences in numerical order to make everything easier to parse. This is literally just a counting loop.##
touch ${prefix}_${gene}_headers
while [[ $x -le $n ]]
do
	if [[ $x -le 9 ]]
 	then 
		echo ">seq_00000${x}" >> ${prefix}_${gene}_headers
		x=$(( $x + 1 ))
	elif [[ $x -le 99 ]] && [[ $x -ge 10 ]]
	then
		echo ">seq_0000${x}" >> ${prefix}_${gene}_headers
		x=$(( $x + 1 ))
	elif [[ $x -le 999 ]] && [[ $x -ge 100 ]]
	then 
		echo ">seq_000${x}" >> ${prefix}_${gene}_headers
		x=$(( $x + 1 ))
	elif [[ $x -le 9999 ]] && [[ $x -ge 1000 ]]
	then
		echo ">seq_00${x}" >> ${prefix}_${gene}_headers
		x=$(( $x + 1 ))
	elif [[ $x -le 999999 ]]
	then
		echo ">seq_0${x}" >> ${prefix}_${gene}_headers
		x=$(( $x + 1 ))
	elif [[ $x -ge 100000 ]]
	then
		echo ">seq_${x}" >> ${prefix}_${gene}_headers
		x=$(( $x + 1 ))
	fi
done
##combine your list of fasta headers with the unique sequences you organized earlier into a fasta-formatted file with every unique ASV from your dataset.##
paste -d '\n' ${prefix}_${gene}_headers temp_seqs > ${prefix}_${gene}_combined_ASVs.fasta
rm temp_seqs ${prefix}_${gene}_headers #clean up interstitial files##
cd ${dirr}

##get or make taxfile. This process looks for an NCBI taxonomy file in your project directory, and if it can't find one, it installs and downloads one.#
cd ${dirr}
if ls ncbi*.csv* 1> /dev/null 2>&1; then
	echo "ncbi tax file found, beginning tax assessment"
	cp ncbi*.csv* ${dirr}/${gene}_out/
	cd ${dirr}/${gene}_out/
	gunzip ncbi*.csv.gz
else
	echo "no ncbi tax file found, attempting to install and run ncbitax2lin."
	pip install -U ncbitax2lin
	wget -N ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
	mkdir -p taxdump && tar zxf taxdump.tar.gz -C ./taxdump
	ncbitax2lin --nodes-file taxdump/nodes.dmp --names-file taxdump/names.dmp
	cp ncbi*.csv* ${dirr}/${gene}_out/
	cd ${dirr}/${gene}_out/
 	gunzip ncbi*.csv.gz
  fi


cd ${dirr}/${gene}_out
ncbi=$( ls ncbi*.csv | head -n1 | awk '{print $1}')
echo "now doing local blast search, this may some time. You have $n sequences to align. You can check your progress in the ${gene}_out folder with the command: tail ${prefix}_${gene}_raw_blast_out"
echo "localdat is $localdat, prefix is $prefix"
##perform local BLAST according to your specifications. You may need to alter the thread_count for your cluster. You can reduce the culling limit to 1 to really limite your blast hits to the top 1 hit.##
if [[ ${return_low} == "TRUE" ]]
then
	echo "you set return_low to TRUE, so BLAST will return the top bitscore matches regardless of percent identity."
  	blastn -db ${dirr}/${db_dirr}/${localdat} -query ${dirr}/${gene}_out/${prefix}_${gene}_combined_ASVs.fasta -outfmt "6 qseqid sacc pident length stitle bitscore" -culling_limit 10 -num_threads 4 -out ${prefix}_${gene}_raw_blast_out
else
 	echo "you set return_low to FALSE, or did not enter a valid TRUE/FALSE value, so BLAST will only return hits above ${cutoff} percent identity, regardless of score."		
  	blastn -db ${dirr}/${db_dirr}/${localdat} -query ${dirr}/${gene}_out/${prefix}_${gene}_combined_ASVs.fasta -outfmt "6 qseqid sacc pident length stitle bitscore" -culling_limit 10 -num_threads 4 -out ${prefix}_${gene}_raw_blast_out -perc_identity ${cutoff}
fi

##define your raw outfiles and remove empty lines which represent no-hits##
blastout=${prefix}_${gene}_raw_blast_out
sed -i '/^#/d' $blastout
echo "blast done, making outfiles"

#make a list of no-hits#
grep ">" ${prefix}_${gene}_combined_ASVs.fasta | awk -F">" '{print $2}' | sort > out1 #a list of every unique sequence, just the headers from your input fasta file
cut -f1 ${blastout}| sort | uniq > out2 #grab the sequence numbers from the blast output, remember no-hits are simply not returned in tab-blast results so this is a list of every sequence with at least one hit#
comm -23 out1 out2 > list_of_no_hits ##sequences that are in the full fasta but not the blast output have no hits in the local blast, they are added to this list##
totalseqs=$( wc -l out1 | awk '{print $1}') ##total number of sequences in your input file##
totalhits=$( wc -l out2 | awk '{print $1}') ##total number of sequencs your local BLAST found a hit for##
nohits=$( wc -l list_of_no_hits | awk '{print $1}') ##number of sequences with no local BLAST hit##
echo "there were ${totalhits} raw BLAST hits out of ${totalseqs} unique sequences, and ${nohits} raw BLAST no-hits. If this number seems too high, consider altering your filtering parameters, changing the blast parameters or adding taxa to your reference database."

##there should always be SOME blast hits. If there aren't any, you probably had a problem with your reference database, usually the program can't find it for some reason. Exit with a fail state.##
if [[ ${totalhits} -eq 0 ]]
then
	echo "there were no local blast hits. This is a problem. Exiting. Check your query fasta, ${prefix}_${gene}_combined_ASVs.fasta and your seqs.txt files for errors."
     	exit 1;
fi
echo "done with local blast, now doing remote blast. There are ${n} sequences to align. This may take many hours. This option is not recommended if you have >50,000 sequences to align. Do NOT set taxa_rra to 0 and then choose this option."
##do the same as above, but for the remote blast. Note: no thread_num because remote.##
if [[ ${return_low} == "TRUE" ]]
then
	echo "you set return_low to TRUE, so BLAST will return the top bitscore matches regardless of percent identity."
  	blastn -db nt -query ${dirr}/${gene}_out/${prefix}_${gene}_combined_ASVs.fasta -outfmt "6 qseqid sacc pident length stitle bitscore staxids" -culling_limit 10 -out ${prefix}_${gene}_remote_raw_blast_out -remote
else
 	echo "you set return_low to FALSE, or did not enter a valid TRUE/FALSE value, so BLAST will only return hits above ${cutoff} percent identity, regardless of score."		
  	blastn -db nt -query ${dirr}/${gene}_out/${prefix}_${gene}_combined_ASVs.fasta -outfmt "6 qseqid sacc pident length stitle bitscore staxids" -culling_limit 10 -out ${prefix}_${gene}_remote_raw_blast_out -perc_identity ${cutoff} -remote
  fi
#make a list of your unique sequences

remote_blastout=${prefix}_${gene}_remote_raw_blast_out
sed -i '/^#/d' $blastout
echo "blast done, making outfiles"

#make a list of no-hits#
cut -f1 ${remote_blastout}| sort | uniq > out2
comm -23 out1 out2 > remote_list_of_no_hits
totalseqs=$( wc -l out1 | awk '{print $1}')
totalhits=$( wc -l out2 | awk '{print $1}')
nohits=$( wc -l remote_list_of_no_hits | awk '{print $1}')
echo "there were ${totalhits} raw BLAST hits out of ${totalseqs} unique sequences, and ${nohits} remote BLAST no-hits. If this number seems too high, consider altering your filtering parameters, changing the blast parameters or adding taxa to your reference database."
cp out1 temp_seqlist
rm out1 out2


##add your species and taxid to your hits using the NCBI taxonomy file##
echo "adding species and taxid to your remote blast results"
cut -f5 $remote_blastout | awk -F" " '{print $1,$2}' > temp_spec ##species, from the NCBI raw blast results. Imperfect, sometimes will get PREDICTED or "severe acute" or something##
cut -f7 $remote_blastout > temp_taxids #taxids are included in the remote blast results#
cut -f6 $remote_blastout > temp_scores #grab bitscores while we're at it#
cut -f1-4 $remote_blastout | paste - temp_spec temp_taxids temp_scores > ${remote_blastout}_with_tax ##all this really does is slightly reorganize the results and simplify the species/description in the raw outfile##
rm temp_spec temp_taxids temp_scores

echo "adding species and taxid to your local blast results"
cut -f5 $blastout | awk -F" " '{print $1,$2}' > temp_spec
cut -f5 $blastout | awk -F"taxid=" '{print $2}' | awk -F" " '{print $1}' > temp_taxids ##we added the taxid to our local sequences in step 1 to make this process easier. This grabs that value from the description. This is the only step where this is different from the remote BLAST##
cut -f6 $blastout > temp_scores
cut -f1-4 $blastout | paste - temp_spec temp_taxids temp_scores > ${blastout}_with_tax
rm temp_spec temp_taxids temp_scores

##modify your ncbi tax file to contain only taxa within your reference database##
echo "now modifying your taxonomy file to limit your search space. This saves time."
cut -f6 ${blastout}_with_tax | sort | uniq > all_local_taxids ##get a list of all the taxids, both local and remote##
cut -f6 ${remote_blastout}_with_tax | sort | uniq > all_remote_taxids
cat all_local_taxids all_remote_taxids | sort | uniq > all_taxids
sed -i 's/$/,/g' all_taxids ##format this file just in case##
sed -i 's/^/\^/g' all_taxids ##remove special characters from this file for formatting reasons##
grep -f all_taxids $ncbi | cut -f1,3,4,5,6,7,8 -d"," > ${ncbi}_r ##extract from the full NCBI taxonomy file only the taxids present in the raw blast outfiles. Because we search this file a lot, this is way faster than using the full database, which contains info for EVERY SPECIES IN NCBI.
sed -i 's/,/\t/g' ${ncbi}_r
rm all_remote_taxids all_local_taxids all_taxids
tot=$( wc -l temp_seqlist | awk '{print $1}')

 echo "making your taxonomic assignment files and beginning to assign jobs."
tot_per_file=$( awk -v a1=$tot -v a2=$max_jobs 'BEGIN { x+=(a1/a2); printf("%.0f", (x == int(x)) ? x : int(x)+1) }' ) ##this looks stupid but it's only because bash hates floating point arithmetic. This literally just divides your max_jobs by the number of sequences we need to assign taxonomy for, and then rounds it appropriately to determine how many jobs to submit##
if [[ ${tot_per_file} -eq 0 ]]
	then
 		 tot_per_file=1
	fi
echo "there were $tot samples to assign taxonomy for and $tot_per_file sample(s) per job."
x=1
##same loop as before, just add sequences (seq_0001 through seq_000n) to lists based on how many sequences/job and how many jobs. Just walks through the list of sequences chunking it into x files.
while [[ $x -le ${max_jobs} ]];
do
	if [[ -s temp_seqlist ]];
  	then
   		head -n ${tot_per_file} temp_seqlist > seqlist_${x}
    		sed -i "1,${tot_per_file}d" temp_seqlist
   		x=$(( $x + 1 ))
 	else
  		x=$(( $max_jobs + 1 ))
  	fi
done
rm temp_seqlist
tot2=$(( $tot * 2 )) ##double this number because we submit two jobs, one remote and one local, for each sequence##
###this loop, like in steps 2-5, checks for an outfile, submits a job, checks that the job is running, waits for all jobs to finish, and then checks to ensure the correct number of outfiles. This is hopefully overkill on many clusters, but necessary on RCDS where job failures are common.##
while [[ $num_outs -ne $tot2 ]];
do
	for fil in seqlist_*;
 	do
  		x=$( echo $fil | awk -F"_" '{print $2}')
    		while true;
     		do
     			if [[ ! -s ${prefix}_${gene}_best_blast_hits.out_${x} ]]; ##check if local outfile exists, or submit##
	 		then
	 			echo "outfile for $fil does not yet exist or is empty. Doing $fil."
     				res=$(sbatch ${dirr}/scripts/run_tax.sh $x $prefix $gene $tot_per_file $blastout $ncbi $dirr $remote_blastout)
   				if [[ ! -s remote_${prefix}_${gene}_best_blast_hits.out_${x} ]];
       				then
       					sbatch ${dirr}/scripts/run_tax_remote.sh $x $prefix $gene $tot_per_file $remote_blastout $ncbi $dirr
       				fi
	   			if squeue -u $user | grep -q "${res##* }"; 
   				then
   					echo "job ${res##* } for $fil submitted successfully."
       					break
       				elif [[ -f tax.${res##* }.err ]];
	  			then
	  				echo "job ${res##* } for $fil submitted successfully."
      					break
      				else
	  				echo "job ${res##* } did not submit. Trying again."
				fi
     			elif [[ ! -s remote_${prefix}_${gene}_best_blast_hits.out_${x} ]]; ##if local exists, double check that remote outfile ALSO exists, or submit##
			then
   				sbatch ${dirr}/scripts/run_tax_remote.sh $x $prefix $gene $tot_per_file $blastout $ncbi $dirr
       				break
			else
	 			echo "all samples from $fil already done."
     				break
      			fi
    		done
       	done	
	while true; ##wait for jobs to run, check when done##
     	do
       		sleep 3s
	 	ck="squeue -u ${user}"
		chck=$($ck)
  		check=$(echo "$chck" | grep "tax" | wc -l | awk '{print $1}')
		echo "waiting for jobs to run. There are $check jobs left"
       		if [[ $check -eq 0 ]];
	 	then
        		echo "no jobs left, checking that jobs ran successfully." 
         		break
       		fi 
	done
  	cat *_best_blast_hits.out_* > outslist
   	num_outs=$( wc -l outslist | awk '{print $1}')
    	if [[ $num_outs -gt $tot2 ]]
     	then
      		num_outs=$tot2
	fi
    	echo "there are $tot2 sequences to assign and $num_outs sequences successfully assigned. If these numbers match, moving on to taxonomy. If not, checking each sample and re-submitting jobs as needed."
done
    
#cat your files and make a header for the best hits table.


cat ${prefix}_${gene}_best_blast_hits.out_* | cut -f1-13 | sort -k1 > ${prefix}_${gene}_best_blast_hits.out
cat remote_${prefix}_${gene}_best_blast_hits.out_* | cut -f2-13 | sort -k1 > remote_${prefix}_${gene}_best_blast_hits.out
eval "$(conda shell.bash hook)"
conda activate $env_name
pref=$CONDA_PREFIX
Rscript ${dirr}/scripts/join.R ${pref}/lib/R/library ${prefix}_${gene}_best_blast_hits.out remote_${prefix}_${gene}_best_blast_hits.out ${prefix} ${gene} $PWD ##literally a script just to do inner_join() in case the remote and local outfiles aren't in the same order, possible if jobs failed halfway or something##

echo "done with choosing best blast hits, now creating and formatting outfiles."
	
#echo "sequence	seqnum	local_identity	local_species	local_taxid	local_phylum	local_class	local_order	local_family	local_genus	local_bitscore	local_num_spec_in_best_hit	local_all_spec_in_best_hit	remote_identity	remote_species	remote_taxid	remote_phylum	remote_class	remote_order	remote_family	remote_genus	remote_bitscore	remote_num_spec_in_best_hit	remote_all_spec_in_best_hit" > ${prefix}_${gene}_best_blast_hits.header
#mv ${prefix}_${gene}_best_blast_hits.out2 > ${prefix}_${gene}_best_blast_hits.txt

#clean up outfiles
rm *${prefix}_${gene}_best_blast_hits.out*
rm seqlist_*
rm tax.*.err
rm tax.*.out
