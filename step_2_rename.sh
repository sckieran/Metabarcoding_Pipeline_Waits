#!/bin/bash

nin=your_barcodes.txt

for fil in *_R1.fastq
do
	base=$(echo $fil | awk -F"_" '{print $1}')
	if grep -q "$base" $nin
	then
	sample=$( grep -m1 "$base" $nin | awk '{print $1}')
	echo "$fil $base $sample"
	mv $fil ${sample}_R1.fastq
	mv ${base}_R2.fastq ${sample}_R2.fastq
	fi
done
