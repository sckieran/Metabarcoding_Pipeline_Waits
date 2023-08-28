#!/bin/bash

#SBATCH -J ttb_j
#SBATCH -e ttb.%j.err
#SBATCH -o ttb.%j.out


infil=$1
dir=$2
gene=$3
prefix=$4

y=$( echo $infil | awk -F"_" '{print $2}')

cd ${dir}/${gene}_out

touch ${prefix}_${gene}_taxatable.txt_${y}

while read fil;
do
  base=$(echo $fil | awk -F"_" '{print $1}')
	x=1
	n=$( wc -l $fil | awk '{print $1}')
	while [[ $x -le $n ]]
	do
			ln="sed -n ${x}p $fil"
			lin=$($ln)
			seq=$(echo $lin | awk  '{print $1}')
			reads=$(echo $lin | awk  '{print $2}')
			taxline=$(grep -w  -m1 "$seq" ${prefix}_${gene}_best_blast_hits.txt | awk -v OFS='\t' '{print $3,$4" "$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16" "$17,$18,$19,$20,$21,$22,$23}')
   			sp2=$(grep -w  -m1 "$seq" ${prefix}_${gene}_best_blast_hits.txt | cut -f13)
      			if cat ${prefix}_${gene}_taxatable.txt_${y} | grep -q -w "$base	$seq";
	 		then
    				echo "this seq already in taxatable for $base."
			else
				echo "$base	$seq	$reads	$taxline	$sp2" >> ${prefix}_${gene}_taxatable.txt_${y}
    			fi
			x=$(( $x + 1 ))
	done
done < $infil
