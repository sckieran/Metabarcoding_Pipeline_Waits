#!/bin/bash
#

retmax=10
comb="sep"

while getopts ":n:t:g:d:r:c:" opt; do
  case $opt in
    n) prefix="$OPTARG"
    ;;
    t) taxlist="$OPTARG"
    ;;
    g) genelist="$OPTARG"
    ;;
    d) dirr="$OPTARG"
    ;;
    r) retmax="$OPTARG"
    ;;
    c) comb="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    exit 1
    ;;
  esac
done
#make the output directory, if it doesn't exist#
mkdir -p $dirr

module load R/4.2.3
module load ncbi-blast
#run rentrez#
Rscript ${dirr}/scripts/query_rentrez.R -n $prefix -t $taxlist -g $genelist -d $dirr -r $retmax 

echo "done downloading reference sequences. Moving on to add taxids to reference sequences."
cp $genelist $dirr
cp $taxlist $dirr
cd $dirr

head -n1 $genelist | sed "s:\t:\n:g" > list_of_genes.txt

for fil in *_sequences.fasta
do
	grep ">" $fil > temp_lines
	taxid=$(head -n1 $fil | awk -F"taxid=" '{print $2}')
	echo "adding taxid $taxid to $fil"
	while read p;
	do
		sed -i "s@${p}@${p}\ttaxid=${taxid}@g" $fil
	done <temp_lines
done
echo "done processing fastas, building database."
if [ $comb == "comb" ]
then
	cat *_sequences.fasta > ${prefix}_all_genes_sequence_database.fasta
	makeblastdb -dbtype nucl -in ${prefix}_all_genes_sequence_database.fasta -out ${prefix}_reference  -parse_seqids -blastdb_version 5
	mkdir -p "all_genes_ref_fastas"
	mv *_sequences.fasta ./all_genes_ref_fastas/
else
	while read p; 
	do
  	mkdir -p "$p"
	cat *_${p}_sequences.fasta > ${prefix}_${p}_sequence_database.fasta
	makeblastdb -dbtype nucl -in ${prefix}_${p}_sequence_database.fasta -out ${prefix}_${p}_reference  -parse_seqids -blastdb_version 5
	mv *_${p}_sequences.fasta ./${p}
	done < list_of_genes.txt
fi


