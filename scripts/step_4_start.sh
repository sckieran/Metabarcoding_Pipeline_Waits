#!/bin/bash
#

pattern1="_R1.fastq"
pattern2="_R2.fastq"
params1="params_file"
score="bitscore"

while getopts ":n:g:d:m:p:q:r:l:b:s:" opt; do
  case $opt in
    n) prefix="$OPTARG"
    ;;
    g) genelist="$OPTARG"
    ;;
    d) dirr="$OPTARG"
    ;;
    m) params1="$OPTARG"
    ;;
    p) pattern1="$OPTARG"
    ;;
    q) pattern2="$OPTARG"
    ;;
    r) db_dirr="$OPTARG"
    ;;
    l) rlib="$OPTARG"
    ;;
    b) localdat="$OPTARG"
    ;;
    c) cutoff="$OPTARG"
    ;;
    t) return_low="$OPTARG"
    ;;
    s) score="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    exit 1
    ;;
    esac
done

module load R/4.2.3
module load ncbi-blast


if [[ "$score" = "bitscore" ]] && [[ "$local" = "remote" ]]
then
echo "bitscore-based assessment not yet available for remote BLAST. Please select local BLAST or pident-based assessment."
fi


cd ${dirr}
head -n1 $genelist | sed 's/\t/\n/g' > list_of_genes.txt

while read p;
do
	params=${params1}_${p}
	gene=${p}
	if [[ -z ${localdat} ]]
	then
		localdat=${prefix}_${gene}_reference
	fi
 	if [[ -z ${db_dirr} ]]
	then
		db_dirr=${dirr}/reference_database
	fi
	mkdir -p ${gene}_dada_out ${gene}_dada_out/filtered/
	echo "gene is $gene"
 	project_dir=$(sed -n 1p $params | awk -F"=" '{print $2}' | awk -F" #" '{print $1}')
	rra_cutoff=$(sed -n 2p $params | awk -F"=" '{print $2}' | awk -F" #" '{print $1}') #relative read abundance cutoff. Default is 0.001 (0.1%). That means if a SAMPLE has 10,000 reads total, this removes any ASV with <10 reads. It is reasonably strict if you have high read counts.
	filt_dir=$(sed -n 3p $params | awk -F"=" '{print $2}' | awk -F" #" '{print $1}') #path to filtered fastas will be /path/to/your/dir/your_gene/filtered/, path to outfiles will be /path/to/your/dir/your_gene_dada_out/
	multi=$(sed -n 11p $params | awk -F"=" '{print $2}' | awk -F" #" '{print $1}') #true or false for multithreading in DADA2

 	#build your filtering parameters based on your params_file
 	bash ${dirr}/scripts/step_3_p1_make_filt_parameters.sh ${params}
	cat ${dirr}/scripts/filter_dada_part1.R filtline ${dirr}/scripts/filter_dada_part2.R > ${dirr}/scripts/filter_and_process_dada.R
	rm filtcomm filtline
	
 	##this script will randomly subsample 50 of your samples down to the first 100k reads to speed up the error learning process. This is because the Rscript was throwing an unusual error that I think is related to the very large processing needs of learnerrors when it uses very large files.
	mkdir -p ${dirr}/${gene}_dada_out/subsampled ${dirr}/results_tables
	#this script will take a params file that contains your gene name and filtering parameters you provide it to process your reads through dada2.
	echo "starting the dada process for ${gene}"
	
 	Rscript ${dirr}/scripts/filter_and_process_dada.R -z ${rlib} -d ${project_dir} -g ${gene} -p ${pattern1} -q ${pattern2} -n ${prefix} -e ${rra_cutoff} -m ${multi}

	cd ${gene}_dada_out
	echo "done with DADA2 processing. Beginning blast and taxtable."
	cat *_seqs.txt | cut -f1 -d"," | sort | uniq > temp_seqs
	sed -i '/^$/d' temp_seqs

	##make query fasta from seqlist#
	x=1
	n=$(wc -l temp_seqs | awk '{print $1}')
	touch ${prefix}_${gene}_combined_ASVs.fasta
	while [[ $x -le $n ]]
	do
		in="sed -n ${x}p temp_seqs"
		seq=$($in)
		if [[ $x -le 9 ]]
		then 
			echo ">seq_0000${x}" >> ${prefix}_${gene}_combined_ASVs.fasta
			echo "$seq" >> ${prefix}_${gene}_combined_ASVs.fasta
			x=$(( $x + 1 ))
		elif [[ $x -le 99 ]] && [[ $x -ge 10 ]]
		then
			echo ">seq_000${x}" >> ${prefix}_${gene}_combined_ASVs.fasta
			echo "$seq" >> ${prefix}_${gene}_combined_ASVs.fasta
			x=$(( $x + 1 ))
		elif [[ $x -le 999 ]] && [[ $x -ge 100 ]]
		then 
			echo ">seq_00${x}" >> ${prefix}_${gene}_combined_ASVs.fasta
 			echo "$seq" >> ${prefix}_${gene}_combined_ASVs.fasta
			x=$(( $x + 1 ))
		elif [[ $x -le 9999 ]] && [[ $x -ge 1000 ]]
		then
			echo ">seq_0${x}" >> ${prefix}_${gene}_combined_ASVs.fasta
			echo "$seq" >> ${prefix}_${gene}_combined_ASVs.fasta
			x=$(( $x + 1 ))
		elif [[ $x -ge 10000 ]]
		then
			echo ">seq_${x}" >> ${prefix}_${gene}_combined_ASVs.fasta
			echo "$seq" >> ${prefix}_${gene}_combined_ASVs.fasta
			x=$(( $x + 1 ))
		fi

	done
	rm temp_seqs
	cd ${dirr}

	##get or make taxfile#
	cd ${dirr}
	if ls ncbi*.csv* 1> /dev/null 2>&1; then
		echo "ncbi tax file found, beginning tax assessment"
		cp ncbi*.csv* ${dirr}/${gene}_dada_out/
		cd ${dirr}/${gene}_dada_out/
		gunzip ncbi*.csv.gz
	else
 		echo "no ncbi tax file found, attempting to install and run ncbitax2lin."
		pip install -U ncbitax2lin
		wget -N ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
		mkdir -p taxdump && tar zxf taxdump.tar.gz -C ./taxdump
		ncbitax2lin --nodes-file taxdump/nodes.dmp --names-file taxdump/names.dmp
		cp ncbi*.csv* ${dirr}/${gene}_dada_out/
 		cd ${dirr}/${gene}_dada_out/
       		gunzip ncbi*.csv.gz
	fi


if [[ "$score" = "bitscore" ]] && [[ "$local" = "local" ]]
then

bash ${dirr}/scripts/step_4_by_score_local.sh -n ${prefix} -g ${gene} -d ${dirr} -r ${db_dirr} -b ${localdat} -c ${cutoff} -t ${return_low} -s ${score}

elif [[ "$score" = "pident" ]] && [[ "$local" = "local" ]]
then
bash ${dirr}/scripts/step_4_by_pident_local.sh -n ${prefix} -g ${gene} -d ${dirr} -r ${db_dirr} -b ${localdat} -c ${cutoff} -t ${return_low} -s ${score}

elif [[ "${score}" = "pident" ]] && [[ "$local" = "remote" ]]
then
bash ${dirr}/scripts/step_4_by_pident_remote.sh -n ${prefix} -g ${gene} -d ${dirr} -r ${db_dirr} -b ${localdat} -c ${cutoff} -t ${return_low} -s ${score}
fi

done < list_of_genes.txt
