#!/bin/bash
#

conda activate metab_pipeline_sblair

while getopts ":n:t:g:d:r:h:l:k:" opt; do
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
    h) db_dirr="$OPTARG"
    ;;
    l) rlib="$OPTARG"
    ;;
    k) key="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    exit 1
    ;;
  esac
done
#make the output directory, if it doesn't exist#
if [[ -z ${retmax} ]]
then
retmax=10
fi
if [[ -z ${db_dirr} ]]
then
db_dirr=reference_database
fi
mkdir -p $db_dirr

#run rentrez#
Rscript ${dirr}/scripts/query_rentrez.R -z ${rlib} -n $prefix -t $taxlist -g $genelist -d ${dirr}/$db_dirr -r $retmax -k $key

echo "done downloading reference sequences. If you have additional off-target sequences to add here, simply ensure that they are in this folder, that each file name ends *_gene1_sequences.fasta and that each sequence line begins with a '>' and a valid NCBI accession number to continue. You must have a different copy of these fastas for each gene1...geneN that you are building a database for, the program will only include fastas for a gene's reference if it ends *_gene1_sequences.fasta."
