#!/bin/bash
#
while getopts ":n:g:d:c:l:" opt; do
  case $opt in
    n) prefix="$OPTARG"
    ;;
    g) genelist="$OPTARG"
    ;;
    d) dirr="$OPTARG"
    ;;
    c) cutoff="$OPTARG"
    ;;
    l) rlib="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    exit 1
    ;;
  esac
done


module load R/4.2.3
module load ncbi-blast
cd ${dirr}

head -n1 $genelist | sed 's/\t/\n/g' > list_of_genes.txt
while read p;
do
	gene=${p}
	cd ${dirr}/${gene}_dada_out

	##get or make taxfile#
	echo "using your identity percent cutoff ${cutoff} to re-run low-identity sequences."
	awk -v c=$cutoff '($3 < c)' ${prefix}_${gene}_best_blast_hits.txt | cut -f1 > rerun_list
	grep -A1 -f rerun_list ${prefix}_${gene}_combined_ASVs.fasta > ${prefix}_${gene}_reruns.fasta
	sed -i 's/--//g' ${prefix}_${gene}_reruns.fasta
	sed -i '/^$/d' ${prefix}_${gene}_reruns.fasta
	rm rerun_list
	num_reruns=$( grep ">" ${prefix}_${gene}_reruns.fasta | wc -l | awk '{print $1}')

	echo "blasting ${num_reruns} sequences using remote NCBI database, this may take some time."
	ncbi=$( ls ncbi*.csv | head -n1 | awk '{print $1}')
	blastn -query  ${prefix}_${gene}_reruns.fasta -db nt -out ${prefix}_${gene}_remote_blast_out -outfmt '6 qseqid sacc pident length stitle staxids' -remote -max_target_seqs 20
	blastout=${prefix}_${gene}_remote_blast_out
	echo "blast done, making outfiles"

	#make your outfile
	echo "seqnum	identity	species	taxid	phylum	class	order	family	genus" > ${prefix}_${gene}_best_remote_blast_hits.out

	#make a list of no-hits#
	grep ">" boco_12S_reruns.fasta | awk -F">" '{print $2}' | sort  > out1
	cut -f1 ${prefix}_${gene}_remote_blast_out | sort | uniq > out2
	comm -23 out1 out2 > list_of_no_hits
	totalseqs=$( wc -l out1 | awk '{print $1}')
	totalhits=$( wc -l out2 | awk '{print $1}')
	nohits=$( wc -l list_of_no_hits | awk '{print $1}')
	echo "there were ${totalhits} hits out of ${totalseqs} unique sequences, and ${nohits} no-hits. If this number seems too high, consider altering your filtering parameters, changing the blast parameters or adding taxa to your reference database."

	#make a list of your unique sequences
	cp out1 temp_seqlist
	rm out1 out2

	##modify your ncbi tax file to contain only taxa within your reference database##
	echo "now modifying your taxonomy file to limit your search space. This saves time."
	cut -f6 ${blastout} | sort | uniq > all_taxids
	sed -i 's/$/,/g' all_taxids
	sed -i 's/^/\^/g' all_taxids
	grep -f all_taxids $ncbi | cut -f1,3,4,5,6,7,8 -d"," > ${ncbi}_r
	sed -i 's/,/\t/g' ${ncbi}_r
	tot=$( wc -l temp_seqlist | awk '{print $1}')
 
	echo "adding species and taxid to your blast results"
	cut -f6 ${blastout} > ${blastout}_edit
	cut -f6 ${blastout} | sort | uniq > all_taxids
	while read t;
	do
		t=$( echo ${t} | awk -F";" '{print $1}')
		spec=$( grep -w "$t" ${ncbi}_r | awk -F"\t" '{print $7}')
		sed -i "s/${t}$/${t}\t${spec}/g" ${blastout}_edit
	done < all_taxids
	cut -f1-4 $blastout | paste - ${blastout}_edit > ${blastout}_with_tax
	rm ${blastout}_edit all_taxids

	##go through your list of sequence hits one by one
	echo "now evaluating your best hits and adding detailed taxonomy information"
	while read p;
	do
		if cat list_of_no_hits | grep -w -q "${p}"
		then
			echo "no hit found by BLAST for ${p}".
			echo "${p}	0	No Hit	NA	NA	NA	NA	NA	NA" >> ${prefix}_${gene}_best_remote_blast_hits.out
		else
			num=$(grep -n "^${p}$" temp_seqlist | awk -F":" '{print $1}')
			echo "doing ${num} (${p}) of ${tot}"
			grep -w "${p}" ${blastout}_with_tax | sort -k3 -nr > temp_seq #get all the identities for the hits for a single sequence
			n=$(cut -f3 temp_seq| uniq -c | head -n1 | awk '{print $1}') ##how many of the top identity% are there? If one, pull that and call the hit. If more than one, then go to next.
			if [[ $n -gt 1 ]]
			then
				echo "more than one equally-good blast hit available. Walking up the taxonomy tree to reach consensus."
				id=$(cut -f3 temp_seq | head -n1 | awk '{print $1}')
				grep "$id" temp_seq > temp_choose ##pull all the matches with the top identity # for a sequence into a new file##
				spec_number=$(cut -f6 temp_choose | sort | uniq | wc -l | awk '{print $1}') ##how many species are in the equally good hits? If one, pull that and call the hit. if more than one, then go to next.##
				if [[ $spec_number -eq 1 ]]
				then
					spec=$(head -n1 temp_choose |cut -f6 | awk '{print $1,$2}') 
					echo "one species amongst best hits. species is ${spec}."
					st=$(head -n1 temp_choose | awk -F"\t" -v OFS='\t' '{print $1,$3,$6,$5}')
			 		taxid=$(head -n1 temp_choose | awk -F"\t" '{print $5}')
					tax=$(grep -w "^${taxid}" ${ncbi}_r | awk -v OFS='\t' '{print $2,$3,$4,$5,$6}')	
					echo "${st}	${tax}" >> ${prefix}_${gene}_best_remote_blast_hits.out
				else #if multiple species are present in best hits, go through each one and add the taxonomy. This is slow, and is like this because adding a taxmap to the ref database simply does not seem to work.#
					while read z; do taxid=$(echo ${z} | awk '{print $5}'); grep -w "^${taxid}" ${ncbi}_r | cut -f2-6; done < temp_choose > temp_ids
					paste temp_choose temp_ids > temp_tax
					cp temp_ids test_ids
					rm temp_choose temp_ids
					gen_number=$( cut -f11 temp_tax | sort | uniq | wc -l)
					fam_number=$( cut -f10 temp_tax | sort | uniq | wc -l)
					order_number=$( cut -f9 temp_tax | sort | uniq | wc -l)
					class_number=$( cut -f8 temp_tax | sort | uniq | wc -l)
					phylum_number=$( cut -f7 temp_tax | sort | uniq | wc -l)
				fi
				if [[ ${spec_number} -gt 1 ]] && [[ $gen_number -eq 1 ]]
				then
					gen=$( cut -f11 temp_tax | sort | uniq | awk '{print $1}')
					echo "one genus and ${spec_number} species amongst best hits. taxa is ${gen} sp."
					ln="grep -m1 "$gen" temp_tax"
					lin=$($ln)
					st=$(echo $lin | awk -v OFS='\t' '{print $1,$3}')
					tax=$(echo $lin | awk -v OFS='\t' '{print $8,$9,$10,$11,$12}')
					echo "${st}	${gen} sp.	"NA"	${tax}" >> ${prefix}_${gene}_best_remote_blast_hits.out
				elif [[ ${spec_number} -gt 1 ]] && [[ $fam_number -eq 1 ]]
				then
					fam=$( cut -f10 temp_tax | sort | uniq | awk '{print $1}')
					echo "one family, ${gen_number} genera and ${spec_number} species amongst best hits. taxa is ${fam} sp."
					ln="grep -m1 "$fam" temp_tax"
					lin=$($ln)
					st=$(echo $lin | awk -v OFS='\t' '{print $1,$3}')
					tax=$(echo $lin | awk -v OFS='\t' '{print $8,$9,$10,$11,"NA"}')
					echo "${st}	${fam} sp.	"NA"	${tax}" >> ${prefix}_${gene}_best_remote_blast_hits.out
				elif [[ ${spec_number} -gt 1 ]] && [[ $order_number -eq 1 ]]
				then
					ord=$( cut -f9 temp_tax | sort | uniq | awk '{print $1}')
					echo "one order, ${fam_number} families, ${gen_number} genera and ${spec_number} species amongst best hits. taxa is ${ord} sp."
					st=$(grep -m1 "$ord" temp_tax | awk -v OFS='\t' '{print $1,$3}')
					tax=$(grep -m1 "$ord" temp_tax | awk -v OFS='\t' '{print $8,$9,$10,"NA","NA"}')
					echo "${st}	${ord} sp.	"NA"	${tax}" >> ${prefix}_${gene}_best_remote_blast_hits.out
				elif [[ ${spec_number} -gt 1 ]] && [[ $class_number -eq 1 ]]
				then
					class=$( cut -f8 temp_tax | sort | uniq | awk '{print $1}')
					echo "one class, ${order_number} orders, ${fam_number} families, ${gen_number} genera and ${spec_number} species amongst best hits. taxa is ${class} sp."
					ln="grep -m1 "$class" temp_tax"
					lin=$($ln)
					st=$(echo $lin | awk -v OFS='\t' '{print $1,$3}')
					tax=$(echo $lin | awk -v OFS='\t' '{print $8,$9,"NA","NA","NA"}')
					echo "${st}	${class} sp.	"NA"	${tax}" >> ${prefix}_${gene}_best_remote_blast_hits.out
				elif [[ ${spec_number} -gt 1 ]] && [[ $phylum_number -eq 1 ]]
				then
					phylum=$( cut -f7 temp_tax | sort | uniq | awk '{print $1}')
					echo "one phylum, ${class_number} classes, ${order_number} orders, ${fam_number} families, ${gen_number} genera and ${spec_number} species amongst best hits. taxa is ${phylum} sp."
					ln="grep -m1 "$phylum" temp_tax"
					lin=$($ln)
					st=$(echo $lin | awk -v OFS='\t' '{print $1,$3}')
					tax=$(echo $lin | awk -v OFS='\t' '{print $8,"NA","NA","NA","NA"}')
					echo "${st}	${phylum} sp.	"NA"	${tax}" >> ${prefix}_${gene}_best_remote_blast_hits.out
				elif [[ ${spec_number} -gt 1 ]] && [[ $phylum_number -gt 1 ]]
				then
					echo "multiple phyla present in equally-good BLAST hits. Designating as no-hit."
					st=$(head -n1 temp_tax | awk -v OFS='\t' '{print $1,$3,"No Hit","NA","NA","NA","NA","NA"}')
					echo "${st}" >> ${prefix}_${gene}_best_remote_blast_hits.out
				fi
			else 
				st=$(head -n1 temp_seq | awk -F"\t" -v OFS='\t' '{print $1,$3,$6,$5}')
				taxid=$(head -n1 temp_seq | awk -F"\t" '{print $5}')
				tax=$(grep -w "^${taxid}" ${ncbi}_r | awk -v OFS='\t' '{print $2,$3,$4,$5,$6}')
				echo "${st}	${tax}" >> ${prefix}_${gene}_best_remote_blast_hits.out
				rm temp_seq
			fi
		fi
	done < temp_seqlist
	rm temp_seqlist

	echo "done with choosing best blast hits, now creating and formatting outfiles."
	h1="head -n1 ${prefix}_${gene}_best_remote_blast_hits.out"
	head1=$($h1)
	tail -n +2 ${prefix}_${gene}_best_remote_blast_hits.out > ${prefix}_${gene}_best_remote_blast_hits.out.r
	echo "sequence	${head1}" > ${prefix}_${gene}_best_remote_blast_hits.txt
	while read c;
	do
		seqnum=$( echo ${c} | awk '{print $1}') 
		seq=$( grep -w -A1 ">${seqnum}" boco_12S_reruns.fasta | tail -n1 | awk '{print $1}')
		echo "${seq}	${c}" >> ${prefix}_${gene}_best_remote_blast_hits.txt
	done < ${prefix}_${gene}_best_remote_blast_hits.out.r
	rm ${prefix}_${gene}_best_remote_blast_hits.out ${prefix}_${gene}_best_remote_blast_hits.out.r

	echo "making remote taxatable"
	echo "sample	sequence	reads	identity	species	taxid	phylum	class	order	family	genus" > ${prefix}_${gene}_remote_taxatable.txt
	for fil in ./sample_seqfiles/*_seqs.txt
	do
		echo "making taxtable, doing $fil"
		base=$(echo $fil | awk -F"_F_filt.fastq_seqs.txt" '{print $1}')
		x=2
		n=$( wc -l $fil | awk '{print $1}')
		while [[ $x -le $n ]]
		do
 			ln="sed -n ${x}p $fil"
			lin=$($ln)
			seq=$(echo $lin | awk -F"," '{print $1}')
			reads=$(echo $lin | awk -F"," '{print $2}')
			taxline=$(grep -w "$seq" ${prefix}_${gene}_best_remote_blast_hits.txt | awk -v OFS='\t' '{print $3,$4" "$5,$6,$7,$8,$9,$10,$11}')
			echo "$base	$seq	$reads	$taxline" >> ${prefix}_${gene}_remote_taxatable.txt
			x=$(( $x + 1 ))
		done
	done
	cp ${prefix}_${gene}_remote_taxatable.txt ${dirr}/results_tables/
	rm temp*

 	echo "now comparing remote and local BLAST hits"
 	Rscript ${dirr}/scripts/compare_taxlists.R ${rlib} ${dirr}/${gene}_dada_out ${prefix}_${gene}_best_blast_hits.txt ${prefix}_${gene}_best_remote_blast_hits.txt ${prefix} ${gene}
	cp *_comparative_BLAST_hits.txt ${dirr}/results_tables/
	cd ${dirr}
 	echo "done with ${gene}. You can find comparative BLAST hit tables in your project directory/results_tables."

done < list_of_genes.txt
