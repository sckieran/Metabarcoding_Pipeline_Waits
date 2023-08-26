#!/bin/bash

while getopts ":n:h:g:d:e:f:" opt; do
  case $opt in
    n) prefix="$OPTARG"
    ;;
    h) db_dirr="$OPTARG"
    ;;
    g) genelist="$OPTARG"
    ;;
    d) dirr="$OPTARG"
    ;;
    e) extra="$OPTARG"
    ;;
    f) env_name="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    exit 1
    ;;
  esac
done

if [[ -z db_dirr ]]
then
	db_dirr=reference_database
fi

module load ncbi-blast
mkdir -p $db_dirr
cp $genelist $db_dirr
if [[ ! -z ${extra} ]]
then
	cp ${extra}*.fasta $db_dirr
fi
cd $db_dirr

echo "done processing fastas, building database."
#loop over each gene in your genelist#
head -n1 $genelist | sed "s:\t:\n:g" > list_of_genes.txt
while read p; 
do
        mkdir -p "$p"
        cat *_${p}_sequences.fasta > ${prefix}_${p}_sequence_database.fasta
 	echo "checking for duplicate sequences in database fasta"
	eval "$(conda shell.bash hook)"
	conda activate $env_name
 	python -u ${dirr}/scripts/step_2_p1_rmdups.py  ${prefix}_${p}_sequence_database.fasta $PWD    
 	mv temp_out ${prefix}_${p}_sequence_database.fasta
	makeblastdb -dbtype nucl -in ${prefix}_${p}_sequence_database.fasta -out ${prefix}_${p}_reference  -parse_seqids -blastdb_version 5
        mv *_${p}_sequences.fasta ./${p}
done < list_of_genes.txt

