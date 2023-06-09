#!/bin/bash

while getopts ":n:h:g:d:" opt; do
  case $opt in
    n) prefix="$OPTARG"
    ;;
    h) db_dirr="$OPTARG"
    ;;
    g) genelist="$OPTARG"
    ;;
    d) dirr="$OPTARG"
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
cd $db_dirr

for fil in *_sequences.fasta
do
	grep ">" $fil | awk -F" " '{print $1}' > temp_lines
	echo "adding taxid to $fil"
	sed -i 's/>//g' temp_lines
	while read p; #for each sequence in each fasta, query NCBI's taxid lookup to get taxid, and attach it to the end of the sequence title.
	do
 		taxid=$( curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${p}&rettype=fasta&retmode=xml" | grep TSeq_taxid | cut -d '>' -f 2 | cut -d '<' -f 1 | tr -d "\n" | awk '{print $1}')
  		sed -i "/${p}/ s/$/\ttaxid=${taxid}/" $fil
	done < temp_lines
done
rm temp_lines
echo "done processing fastas, building database."
#loop over each gene in your genelist#
head -n1 $genelist | sed "s:\t:\n:g" > list_of_genes.txt
while read p; 
do
        mkdir -p "$p"
        cat *_${p}_sequences.fasta > ${prefix}_${p}_sequence_database.fasta
 	echo "checking for duplicate sequences in database fasta"
	bash ${dirr}/scripts/step_2_p1_rmdups.sh  ${prefix} ${p}     
	makeblastdb -dbtype nucl -in ${prefix}_${p}_sequence_database.fasta -out ${prefix}_${p}_reference  -parse_seqids -blastdb_version 5
        mv *_${p}_sequences.fasta ./${p}
done < list_of_genes.txt

