#!/bin/bash
#


while getopts ":n:t:g:d:r:h:l:k:p:e:s:" opt; do
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
    p) env_name="$OPTARG"
    ;;
    e) email="$OPTARG"
    ;;
    s) genus_search="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    exit 1
    ;;
  esac
done
cd $dirr
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
cp $taxlist ./${db_dirr}
cp $genelist ./${db_dirr}

#module load R/4.2.3
#module load ncbi-blast
#run rentrez#

eval "$(conda shell.bash hook)"
conda activate $env_name
python -u ${dirr}/scripts/query_rentrez.py $prefix ${PWD}/${db_dirr} $genelist $taxlist $retmax $genus_search $key $email

exit_status=$?
if [ "${exit_status}" -ne 0 ];
then
    echo "exit ${exit_status}"
    exit 1
fi
echo "done downloading reference sequences. If you have additional off-target sequences to add here, simply ensure that they are in this folder, that each file name ends *_gene1_sequences.fasta and that each sequence line begins with a '>' and a valid NCBI accession number to continue. You must have a different copy of these fastas for each gene1...geneN that you are building a database for, the program will only include fastas for a gene's reference if it ends *_gene1_sequences.fasta."
