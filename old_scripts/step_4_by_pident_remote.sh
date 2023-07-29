#!/bin/bash
#
score="pident"

while getopts ":n:g:d:l:c:t:s:" opt; do
  case $opt in
	n) prefix="$OPTARG"
	;;
	g) gene="$OPTARG"
	;;
	d) dirr="$OPTARG"
	;;
	l) rlib="$OPTARG"
	;;
	c) cutoff="$OPTARG"
	;;
	t) return_low="$OPTARG"
	;;
	s) score="$OPTARG"
	;;
	\?) echo "Invalid option -$OPTARG" >&2
	exit 1
	;;
  esac
done

module load R/4.2.3
module load ncbi-blast

##start looping over your genelist##
cd ${dirr}/${gene}_dada_out

#calculate how many sequences you're blasting#
	numseqs=$( grep ">"  ${prefix}_${gene}_combined_ASVs.fasta | wc -l | awk '{print $1}')
	echo "now doing remote blast search on ${numseqs} sequences, this may take several hours."
	blastn -db nt -query ${dirr}/${gene}_dada_out/${prefix}_${gene}_combined_ASVs.fasta -outfmt '6 qseqid sacc pident length stitle staxid' -out ${prefix}_${gene}_remote_blast_out -remote -max_target_seqs 20

	blastout=${prefix}_${gene}_remote_blast_out
	echo "blast done, making outfiles"
	#make your outfile
	echo "seqnum	identity	species	taxid	phylum	class	order	family	genus" > ${prefix}_${gene}_best_remote_blast_hits.out
	###############
	
	#make a list of no-hits#
	grep ">" ${prefix}_${gene}_combined_ASVs.fasta | awk -F">" '{print $2}' | sort  > out1
	cut -f1 ${prefix}_${gene}_remote_blast_out | sort | uniq > out2
	comm -23 out1 out2 > list_of_no_hits
	totalseqs=$( wc -l out1 | awk '{print $1}')
	totalhits=$( wc -l out2 | awk '{print $1}')
	nohits=$( wc -l list_of_no_hits | awk '{print $1}')
	echo "there were ${totalhits} hits out of ${totalseqs} unique sequences, and ${nohits} no-hits. If this number seems too high, consider altering your filtering parameters, changing the blast parameters or adding taxa to your reference database."
	#make a list of your unique sequences
	cp out1 temp_seqlist
	rm out1 out2
	##add your species and taxid to your hits##
	cd ${dirr}
	if ls ncbi*.csv* 1> /dev/null 2>&1; then
		echo "ncbi tax file found, beginning tax assessment"
		cp ncbi*.csv* ${dirr}/${gene}_dada_out/
		cd ${dirr}/${gene}_dada_out/	
		gunzip ncbi*.csv.gz
	else
		echo "no ncbi tax file found, attempting to install and run ncbitax2lin."
		pip install -U ncbitax2lin
		wget -N ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
		mkdir -p taxdump && tar zxf taxdump.tar.gz -C ./taxdump
		ncbitax2lin --nodes-file taxdump/nodes.dmp --names-file taxdump/names.dmp
		cp ncbi*.csv* ${dirr}/${gene}_dada_out/
 		cd ${dirr}/${gene}_dada_out/
		gunzip ncbi*.csv.gz
	fi


	##modify your ncbi tax file to contain only taxa within your reference database##
	echo "now modifying your taxonomy file to limit your search space. This saves time."
	cut -f6 ${blastout} | sort | uniq > all_taxids
	sed -i 's/$/,/g' all_taxids
	sed -i 's/^/\^/g' all_taxids
	grep -f all_taxids $ncbi | cut -f1,3,4,5,6,7,8 -d"," > ${ncbi}_r
	sed -i 's/,/\t/g' ${ncbi}_r
	rm all_taxids
	tot=$( wc -l temp_seqlist | awk '{print $1}')


	echo "adding species and taxid to your blast results"
	cut -f6 ${blastout} | sort | uniq > all_taxids
	cut -f6 ${blastout} > ${blastout}_edit
	while read t;
	do
		t=$( echo ${t} | awk -F";" '{print $1}')
		spec=$( grep -w "$t" ${ncbi}_r | awk -F"\t" '{print $7}')
		sed -i "s/${t}$/${t}\t${spec}" ${blastout}_edit
	done < all_taxids
	cut -f1-4 $blastout | paste - ${blastout}_edit > ${blastout}_with_tax
	rm ${blastout}_edit all_taxids


	###################
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
			top=$(cut -f3 temp_seq | uniq | head -n1 | awk -F"." '{print $1}')
			tp=$(cut -f3 temp_seq | uniq | head -n1 | awk '{print $1}')
			if [[ "$return_low" = "FALSE" ]] && [[ ${top} -lt ${cutoff} ]]
			then
				"the top percent identity for this match is below your threshold ${cutoff}. Because you set return_low to FALSE, calling a no-hit for this ASV, but reporting highest % identity."			echo "${p}      ${tp}       No Hit  NA      NA      NA      NA      NA      NA" >> ${prefix}_${gene}_best_remote_blast_hits.out
			elif [[ $n -gt 1 ]] 
			then
				echo "more than one equally-good blast hit available. Walking up the taxonomy tree to reach consensus."
				id=$(cut -f3 temp_seq | head -n1 | awk '{print $1}')
				grep "$id" temp_seq > temp_choose ##pull all the matches with the top identity # for a sequence into a new file##
				spec_number=$(cut -f6 temp_choose | sort | uniq | wc -l | awk '{print $1}') ##how many species are in the equally good hits? If one, pull that and call the hit. if more than one, then go to next.##
				if [[ $spec_number -eq 1 ]]
				then
					spec=$(head -n1 temp_choose | cut -f6 | awk '{print $1,$2}')
					echo "one species amongst best hits. species is ${spec}."
					st=$(head -n1 temp_choose | awk -F"\t"  -v OFS='\t' '{print $1,$3,$6,$5}')
					taxid=$(head -n1 temp_choose | awk -F"\t" '{print $5}')
					tax=$(grep -w "^${taxid}" ${ncbi}_r | awk -v OFS='\t' '{print $2,$3,$4,$5,$6}')
					echo "${st}	${tax}" >> ${prefix}_${gene}_best_remote_blast_hits.out
				else #if multiple species are present in best hits, go through each one and add the taxonomy. This is slow, and is like this because adding a taxmap to the ref database simply does not seem to work.#
					while read z; do taxid=$(echo ${z} | awk '{print $5}'); grep -w "^${taxid}" ${ncbi}_r | cut -f2-6; done < temp_choose > temp_ids
					paste temp_choose temp_ids > temp_tax
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
					ln="grep -m1 "$ord" temp_tax"
					lin=$($ln)
					st=$(echo $lin | awk -v OFS='\t' '{print $1,$3}')
					tax=$(echo $lin | awk -v OFS='\t' '{print $8,$9,$10,"NA","NA"}')
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
	rm temp_sq ${prefix}_${gene}_best_remote_blast_hits.out ${prefix}_${gene}_best_remote_blast_hits.out.r

	echo "making taxatable"
	echo "sample	sequence	reads	identity	species	taxid	phylum	class	order	family	genus" >  ${prefix}_${gene}_taxatable.txt
	for fil in *_seqs.txt
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
			echo "$base	$seq	$reads	$taxline" >> ${prefix}_${gene}_taxatable.txt
			x=$(( $x + 1 ))
		done
	done
	echo "done with ${gene}. You can find taxa tables and raw ASV tables  in your project directory/reports."

	rm temp*
	mkdir -p sample_seqfiles
	mv *_seqs.txt ./sample_seqfiles

	cd ${dirr}
	cp ${dirr}/${gene}_dada_out/*_taxatable.txt ${dirr}/results_tables/

