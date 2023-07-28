#!/bin/bash

prefix=$1
gene=$2
dirr=$3

echo "sample	sequence	reads	identity	taxa	taxid	phylum	class	order	family	genus	bitscore	tax_num	all_species_in_best_hit	seqnum" >  ${prefix}_${gene}_taxatable.txt
	for fil in *_seqs.txt
	do
		echo "making taxtable, doing $fil"
		base=$(echo $fil | awk -F"_" '{print $1}')
		x=2
		n=$( wc -l $fil | awk '{print $1}')
		while [[ $x -le $n ]]
		do
			ln="sed -n ${x}p $fil"
			lin=$($ln)
			seq=$(echo $lin | awk  '{print $1}')
			reads=$(echo $lin | awk  '{print $2}')
			taxline=$(grep -w  -m1 "$seq" ${prefix}_${gene}_best_blast_hits.txt | awk -v OFS='\t' '{print $3,$4" "$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$2}')
			echo "$base	$seq	$reads	$taxline" >> ${prefix}_${gene}_taxatable.txt
			x=$(( $x + 1 ))
		done
	done
	echo "done with ${gene}. You can find taxa tables and raw ASV tables in your project directory/reports."
	rm temp*
	mkdir ${dirr}/${gene}_out/sample_seqfiles
	mv *_seqs.txt ${dirr}/${gene}_out/sample_seqfiles
	cd ${dirr}
	cp ${dirr}/${gene}_out/*_taxatable.txt ${dirr}/results_tables/
