#!/bin/bash

#SBATCH -J tax
#SBATCH -e tax.%j.err
#SBATCH -o tax.%j.out

x=$1
prefix=$2
gene=$3
tot=$4
blastout=$5
ncbi=$6
dir=$7

cd ${dir}/${gene}_out
echo -n "" > ${prefix}_${gene}_best_blast_hits.out_${x}
while read p;
	do
		sqy=$( grep -w -A1 "$p" ${prefix}_${gene}_combined_ASVs.fasta | tail -n1 | awk '{print $1}')
  		if cat list_of_no_hits | grep -w -q "${p}"
		then
			echo "no hit found by BLAST for ${p}".
			echo "${sqy}	${p}	0	No Hit	NA	NA	NA	NA	NA	NA" >> ${prefix}_${gene}_best_blast_hits.out_${x}
		else
			num=$(grep  -n "^${p}$" seqlist_${x} | awk -F":" '{print $1}')
			echo "doing ${num} (${p}) of ${tot}"
			grep -w "${p}" ${blastout}_with_tax | sort -k8 -nr > temp_seq_${x} #get all the identities for the hits for a single sequence
			n=$(cut -f7 temp_seq_${x} | uniq -c | head -n1 | awk '{print $1}') ##how many of the top identity% are there? If one, pull that and call the hit. If more than one, then go to next.
			top_score=$(cut -f7 temp_seq_${x} | head -n1 | awk '{print $1}')
			awk -v d=$top_score '( $8 >= d )'  temp_seq_${x} | awk -v OFS='\t' '{print $1,$2,$3,$4,$5" "$6,$7,$8}'> temp_choose2_${x}
			top_id2=$( sort -k3 -nr temp_choose2_${x} | cut -f3 | head -n1 | awk '{print $1}')
			echo "there are $n hit(s) at best score for sequence ${p}. The top identity % at highest score is $top_id2"
			top_id=$( echo $top_id2 | awk -F"." '{print $1}')
			if [[ $n -eq 1 ]] 
			then	
				ln="head -n1 temp_seq_${x}"
				num_tx=1
				spec=$( head -n1 temp_seq_${x} | cut -f5 | awk '{print $1,$2}')
				lin=$($ln)
				st=$(echo $lin | awk -v OFS='\t' '{print $1,$3,$5" "$6,$7}')
				taxid=$(echo $lin | awk '{print $7}')
				tax=$(grep -w "^${taxid}" ${ncbi}_r | awk -v OFS='\t' '{print $2,$3,$4,$5,$6}')
				echo "${sqy}	${st}	${tax}	${top_score}	${num_tx}	${spec}" >> ${prefix}_${gene}_best_blast_hits.out_${x}
				rm temp_seq_${x}
			elif [[ $n -gt 1 ]] 
			then
				echo "more than one equally-good blast hit available. Walking up the taxonomy tree to reach consensus."
				id=$(cut -f3 temp_choose2_${x} | sort -nr |  head -n1 | awk '{print $1}')
				grep "$id" temp_choose2_${x} > temp_choose_${x} ##pull all the matches with the top identity # for a sequence into a new file##
				spec_number=$(cut -f5 temp_choose_${x} | sort | uniq | wc -l | awk '{print $1}') ##how many species are in the equally good hits? If one, pull that and call the hit. if more than one, then go to next.##
				cut -f5 temp_choose_${x} | sort | uniq > temp_spec_${x}
				sed -zi 's/\n/,/g' temp_spec_${x}
				spec2=$( cat temp_spec_${x} | awk -F"\t" '{print $1}')
				num_tx=${spec_number}
				if [[ $spec_number -eq 1 ]]
				then
					spec=$(head -n1 temp_choose_${x} |cut -f5 | awk '{print $1,$2}') 
					echo "one species amongst best hits. species is ${spec}."
					ln="head -n1 temp_choose_${x}"
					lin=$($ln)
					st=$(echo $lin | awk -v OFS='\t' '{print $1,$3,$5" "$6,$7}')
					taxid=$(echo $lin | awk '{print $7}')
					tax=$(grep -w "^${taxid}" ${ncbi}_r | awk -v OFS='\t' '{print $2,$3,$4,$5,$6}')	
					echo "${sqy}	${st}	${tax}	${top_score}	${num_tx}	${spec2}" >> ${prefix}_${gene}_best_blast_hits.out_${x}
				else #if multiple species are present in best hits, go through each one and add the taxonomy. This is slow, and is like this because adding a taxmap to the ref database simply does not seem to work.#
					while read z;
					do
						taxid=$(echo ${z} | awk '{print $7}')
						grep -w "^${taxid}" ${ncbi}_r | cut -f2-6 >> temp_ids_${x}
					done < temp_choose_${x}
					paste temp_choose_${x} temp_ids_${x} > temp_tax_${x}
					rm temp_choose_${x} temp_ids_${x}
					gen_number=$( cut -f12 temp_tax_${x} | sort | uniq | wc -l)
					fam_number=$( cut -f11 temp_tax_${x} | sort | uniq | wc -l)
					order_number=$( cut -f10 temp_tax_${x} | sort | uniq | wc -l)
					class_number=$( cut -f9 temp_tax_${x} | sort | uniq | wc -l)
					phylum_number=$( cut -f8 temp_tax_${x} | sort | uniq | wc -l)
				fi
				if [[ ${spec_number} -gt 1 ]] && [[ $gen_number -eq 1 ]]
				then
					gen=$( cut -f12 temp_tax_${x} | sort | uniq | awk '{print $1}')
					echo "one genus and ${spec_number} species amongst best hits. taxa is ${gen} sp."
					ln="grep -m1 "$gen" temp_tax_${x}"
					lin=$($ln)
					st=$(echo $lin | awk -v OFS='\t' '{print $1,$3}')
					tax=$(echo $lin | awk -v OFS='\t' '{print $9,$10,$11,$12,$13}')
					echo "${sqy}	${st}	${gen} sp.	"NA"	${tax}	${top_score}	${num_tx}	${spec2}" >> ${prefix}_${gene}_best_blast_hits.out_${x}
				elif [[ ${spec_number} -gt 1 ]] && [[ $fam_number -eq 1 ]]
				then
					fam=$( cut -f11 temp_tax_${x} | sort | uniq | awk '{print $1}')
					echo "one family, ${gen_number} genera and ${spec_number} species amongst best hits. taxa is ${fam} sp."
					ln="grep -m1 "$fam" temp_tax_${x}"
					lin=$($ln)
					st=$(echo $lin | awk -v OFS='\t' '{print $1,$3}')
					tax=$(echo $lin | awk -v OFS='\t' '{print $9,$10,$11,$12,"NA"}')
					echo "${sqy}	${st}	${fam} sp.	"NA"	${tax}	${top_score}	${num_tx}	${spec2}" >> ${prefix}_${gene}_best_blast_hits.out_${x}
				elif [[ ${spec_number} -gt 1 ]] && [[ $order_number -eq 1 ]]
				then
					ord=$( cut -f10 temp_tax_${x} | sort | uniq | awk '{print $1}')
					echo "one order, ${fam_number} families, ${gen_number} genera and ${spec_number} species amongst best hits. taxa is ${ord} sp."
					ln="grep -m1 "$ord" temp_tax_${x}"
					lin=$($ln)
					st=$(echo $lin | awk -v OFS='\t' '{print $1,$3}')
					tax=$(echo $lin | awk -v OFS='\t' '{print $9,$10,$11,"NA","NA"}')
					echo "${sqy}	${st}	${ord} sp.	"NA"	${tax}	${top_score}	${num_tx}	${spec2}" >> ${prefix}_${gene}_best_blast_hits.out_${x}
				elif [[ ${spec_number} -gt 1 ]] && [[ $class_number -eq 1 ]]
				then
					class=$( cut -f9 temp_tax_${x} | sort | uniq | awk '{print $1}')
					echo "one class, ${order_number} orders, ${fam_number} families, ${gen_number} genera and ${spec_number} species amongst best hits. taxa is ${class} sp."
					ln="grep -m1 "$class" temp_tax_${x}"
					lin=$($ln)
					st=$(echo $lin | awk -v OFS='\t' '{print $1,$3}')
					tax=$(echo $lin | awk -v OFS='\t' '{print $9,$10,"NA","NA","NA"}')
					echo "${sqy}	${st}	${class} sp.	"NA"	${tax}	${top_score}	${num_tx}	${spec2}" >> ${prefix}_${gene}_best_blast_hits.out_${x}
				elif [[ ${spec_number} -gt 1 ]] && [[ $phylum_number -eq 1 ]]
				then
					phylum=$( cut -f8 temp_tax_${x} | sort | uniq | awk '{print $1}')
					echo "one phylum, ${class_number} classes, ${order_number} orders, ${fam_number} families, ${gen_number} genera and ${spec_number} species amongst best hits. taxa is ${phylum} sp."
					ln="grep -m1 "$phylum" temp_tax_${x}"
					lin=$($ln)
					st=$(echo $lin | awk -v OFS='\t' '{print $1,$3}')
					tax=$(echo $lin | awk -v OFS='\t' '{print $9,"NA","NA","NA","NA"}')
					echo "${sqy}	${st}	${phylum} sp.	"NA"	${tax}	${top_score}	${num_tx}	${spec2}" >> ${prefix}_${gene}_best_blast_hits.out_${x}
				elif [[ ${spec_number} -gt 1 ]] && [[ $phylum_number -gt 1 ]]
				then
					echo "multiple phyla present in equally-good BLAST hits. Designating as no-hit, but reporting best score and top identity."
					st=$( head -n1 temp_tax_${x} | awk -v OFS='\t' '{print $1,$3,"No Hit","NA","NA","NA","NA","NA"}')
					echo "${sqy}	${st}	${top_score}	${num_tx}	${spec2}" >> ${prefix}_${gene}_best_blast_hits.out_${x}
				fi
		fi	
  	fi
 done < seqlist_${x}
 rm temp_*_${x}
 echo "job seqlist_${x} is done"
 
