#!/bin/bash
#
#
#this script will:
#
#1. look for your data files in the path $dirr/$gene for each gene
#2. produce quality reports.


module load R/4.2.3
num_graphs=24
pattern1="_R1.fastq"
pattern2="_R2.fastq"
##parse command line arguments##
while getopts ":n:g:d:k:p:q:l:" opt; do
  case $opt in
    n) prefix="$OPTARG"
    ;;
    g) gene="$OPTARG"
    ;;
    d) dirr="$OPTARG"
    ;;
    k) num_graphs="$OPTARG"
    ;;
    p) pattern1="$OPTARG"
    ;;
    q) pattern2="$OPTARG"
    ;;
    l) rlib="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    exit 1
    ;;
  esac
done

cd $dirr

##check how many genes you have, although this should be re-written because the behavior is actually identical. Originally I had it so that you didn't need to keep your files in a containing folder called "gene". But it's much easier to organize this way. This just runs the R script for each gene and produces separate plots.
head -n1 $gene | sed "s:\t:\n:g" > list_of_genes.txt

mkdir -p ${dirr}/reports

while read p;
do
	gene=${p}
	mkdir -p ${gene}_dada_out
	Rscript ${dirr}/scripts/quality_check_reads_in_DADA2.R -z ${rlib} -d ${dirr} -g ${gene} -p ${pattern1} -q ${pattern2} -n ${prefix} -k ${num_graphs}
	mv ${dirr}/${gene}/${prefix}*.pdf ${dirr}/reports/
done < list_of_genes.txt



