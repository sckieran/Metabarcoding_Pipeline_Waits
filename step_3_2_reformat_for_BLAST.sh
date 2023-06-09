#!/bin/bash


##makes a list of the unique sequences passing filter in at least one sample##
cat *_seqs.txt | cut -d"," -f1 | sort | uniq | tail -n+2 > your_project_seqs_only.txt

##set input and output names, output will be a fasta file##
in1=your_project_seqs_formatted
inp=your_project_seqs_only.txt


##this loop will create a fasta out of the list of seqs by adding a fasta header line (starts with ">") with seq_1, seq_2, etc. for every seq. Need to do this for BLAST to work.##
x=1
n=$(wc -l $inp | awk '{print $1}')
touch ${in1}.fasta


while [[ $x -le $n ]]
do
in="sed -n ${x}p ${inp}"
seq=$($in)
echo ">seq_${x}" >> ${in1}.fasta
echo "$seq" >> ${in1}.fasta
x=$(( $x + 1 ))
done

