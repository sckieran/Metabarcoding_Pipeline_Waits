#!/bin/bash


out=your_project_taxa_table.txt
lookup=your_project_with_taxa.txt 

echo "sample	sequence	taxa	identity" > $out

for fil in *_seqs.txt
do
	echo "doing $fil"
	base=$(echo $fil | awk -F"_" '{print $1}')
	x=2
	n=$( wc -l $fil | awk '{print $1}')
	while [[ $x -le $n ]]
	do
		ln="sed -n ${x}p $fil"
		lin=$($ln)
		seq=$(echo $lin | awk -F"," '{print $1}')
		reads=$(echo $lin | awk -F"," '{print $2}')
		taxline=$(grep $seq $lookup | awk '{print $1,$2,$3,$4}')
		set -- $taxline
		seq=$1
		taxa="${2} ${3}"
		identity=$4
		echo "$base	$seq	$reads	$taxa	$identity" >> $out
		x=$(( $x + 1 ))
	done
done
		
