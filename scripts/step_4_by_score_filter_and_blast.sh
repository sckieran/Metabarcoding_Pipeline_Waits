#!/bin/bash
#

pattern1="_R1.fastq"
pattern2="_R2.fastq"
params1="params_file"

while getopts ":n:g:d:m:p:q:r:l:b:" opt; do
  case $opt in
    n) prefix="$OPTARG"
    ;;
    g) genelist="$OPTARG"
    ;;
    d) dirr="$OPTARG"
    ;;
    m) params1="$OPTARG"
    ;;
    p) pattern1="$OPTARG"
    ;;
    q) pattern2="$OPTARG"
    ;;
    r) db_dirr="$OPTARG"
    ;;
    l) rlib="$OPTARG"
    ;;
    b) localdat="$OPTARG"
    ;;
    c) cutoff="$OPTARG"
    ;;
    t) return_low="$OPTARG"
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
	params=${params1}_${p}
	gene=${p}
	if [[ -z ${localdat} ]]
	then
		localdat=${prefix}_${gene}_reference
	fi
 	if [[ -z ${db_dirr} ]]
	then
		db_dirr=${dirr}/reference_database
	fi
	mkdir -p ${gene}_dada_out ${gene}_dada_out/filtered/
	echo "gene is $gene"
 	project_dir=$(sed -n 1p $params | awk -F"=" '{print $2}' | awk -F" #" '{print $1}')
	rra_cutoff=$(sed -n 2p $params | awk -F"=" '{print $2}' | awk -F" #" '{print $1}') #relative read abundance cutoff. Default is 0.001 (0.1%). That means if a SAMPLE has 10,000 reads total, this removes any ASV with <10 reads. It is reasonably strict if you have high read counts.
	filt_dir=$(sed -n 3p $params | awk -F"=" '{print $2}' | awk -F" #" '{print $1}') #path to filtered fastas will be /path/to/your/dir/your_gene/filtered/, path to outfiles will be /path/to/your/dir/your_gene_dada_out/
	multi=$(sed -n 11p $params | awk -F"=" '{print $2}' | awk -F" #" '{print $1}') #true or false for multithreading in DADA2

 	#build your filtering parameters based on your params_file
 	bash ${dirr}/scripts/step_3_p1_make_filt_parameters.sh ${params}
	cat ${dirr}/scripts/filter_dada_part1.R filtline ${dirr}/scripts/filter_dada_part2.R > ${dirr}/scripts/filter_and_process_dada.R
	rm filtcomm filtline
	
 	##this script will randomly subsample 50 of your samples down to the first 100k reads to speed up the error learning process. This is because the Rscript was throwing an unusual error that I think is related to the very large processing needs of learnerrors when it uses very large files.
	mkdir -p ${dirr}/${gene}_dada_out/subsampled ${dirr}/results_tables
	#this script will take a params file that contains your gene name and filtering parameters you provide it to process your reads through dada2.
	echo "starting the dada process for ${gene}"
	
 	Rscript ${dirr}/scripts/filter_and_process_dada.R -z ${rlib} -d ${project_dir} -g ${gene} -p ${pattern1} -q ${pattern2} -n ${prefix} -e ${rra_cutoff} -m ${multi}

	cd ${gene}_dada_out
	echo "done with DADA2 processing. Beginning blast and taxtable."
	cat *_seqs.txt | cut -f1 -d"," | sort | uniq > temp_seqs
	sed -i '/^$/d' temp_seqs

	##make query fasta from seqlist#
	x=1
	n=$(wc -l temp_seqs | awk '{print $1}')
	touch ${prefix}_${gene}_combined_ASVs.fasta
	while [[ $x -le $n ]]
	do
		in="sed -n ${x}p temp_seqs"
		seq=$($in)
		if [[ $x -le 9 ]]
		then 
			echo ">seq_0000${x}" >> ${prefix}_${gene}_combined_ASVs.fasta
			echo "$seq" >> ${prefix}_${gene}_combined_ASVs.fasta
			x=$(( $x + 1 ))
		elif [[ $x -le 99 ]] && [[ $x -ge 10 ]]
		then
			echo ">seq_000${x}" >> ${prefix}_${gene}_combined_ASVs.fasta
			echo "$seq" >> ${prefix}_${gene}_combined_ASVs.fasta
			x=$(( $x + 1 ))
		elif [[ $x -le 999 ]] && [[ $x -ge 100 ]]
		then 
			echo ">seq_00${x}" >> ${prefix}_${gene}_combined_ASVs.fasta
 			echo "$seq" >> ${prefix}_${gene}_combined_ASVs.fasta
			x=$(( $x + 1 ))
		elif [[ $x -le 9999 ]] && [[ $x -ge 1000 ]]
		then
			echo ">seq_0${x}" >> ${prefix}_${gene}_combined_ASVs.fasta
			echo "$seq" >> ${prefix}_${gene}_combined_ASVs.fasta
			x=$(( $x + 1 ))
		elif [[ $x -ge 10000 ]]
		then
			echo ">seq_${x}" >> ${prefix}_${gene}_combined_ASVs.fasta
			echo "$seq" >> ${prefix}_${gene}_combined_ASVs.fasta
			x=$(( $x + 1 ))
		fi

	done
	rm temp_seqs
	cd ${dirr}

	##get or make taxfile#
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

	cd ${dirr}/${gene}_dada_out
	ncbi=$( ls ncbi*.csv | head -n1 | awk '{print $1}')
	echo "now doing local blast search"
	blastn -db ${db_dirr}/${localdat} -query ${dirr}/${gene}_dada_out/${prefix}_${gene}_combined_ASVs.fasta -outfmt '6 qseqid sacc pident length stitle bitscore' -out ${prefix}_${gene}_raw_blast_out

	blastout=${prefix}_${gene}_raw_blast_out
	echo "blast done, making outfiles"
	#make your outfile
	echo "seqnum	identity	species	taxid	phylum	class	order	family	genus	top_score	chosen_score" > ${prefix}_${gene}_best_blast_hits.out

	#make a list of no-hits#
	grep ">" ${prefix}_${gene}_combined_ASVs.fasta | awk -F">" '{print $2}' | sort > out1
	cut -f1 ${prefix}_${gene}_raw_blast_out | sort | uniq > out2
	comm -23 out1 out2 > list_of_no_hits
	totalseqs=$( wc -l out1 | awk '{print $1}')
	totalhits=$( wc -l out2 | awk '{print $1}')
	nohits=$( wc -l list_of_no_hits | awk '{print $1}')
	echo "there were ${totalhits} raw BLAST hits out of ${totalseqs} unique sequences, and ${nohits} raw BLAST no-hits. If this number seems too high, consider altering your filtering parameters, changing the blast parameters or adding taxa to your reference database."

	#make a list of your unique sequences
	cp out1 temp_seqlist
	rm out1 out2

	##add your species and taxid to your hits##
	echo "adding species and taxid to your blast results"
	cut -f5 $blastout | awk -F" " '{print $1,$2}' > temp_spec
	cut -f5 $blastout | awk -F"taxid=" '{print $2}' | awk -F" " '{print $1}' > temp_taxids
	cut -f6 $blastout > temp_scores
	cut -f1-4 $blastout | paste - temp_spec temp_taxids temp_scores > ${blastout}_with_tax
	rm temp_spec temp_taxids temp_scores

	##modify your ncbi tax file to contain only taxa within your reference database##
	echo "now modifying your taxonomy file to limit your search space. This saves time."
	cut -f6 ${blastout}_with_tax | sort | uniq > all_taxids
	sed -i 's/$/,/g' all_taxids
	sed -i 's/^/\^/g' all_taxids
	grep -f all_taxids $ncbi | cut -f1,3,4,5,6,7,8 -d"," > ${ncbi}_r
	sed -i 's/,/\t/g' ${ncbi}_r
	rm all_taxids
	tot=$( wc -l temp_seqlist | awk '{print $1}')
	
 	##go through your list of sequence hits one by one
	echo "now evaluating your best hits and adding detailed taxonomy information"
	while read p;
	do
		if cat list_of_no_hits | grep -w -q "${p}"
		then
			echo "no hit found by BLAST for ${p}".
			echo "${p}	0	No Hit	NA	NA	NA	NA	NA	NA" >> ${prefix}_${gene}_best_blast_hits.out
		else
			num=$(grep  -n "^${p}$" temp_seqlist | awk -F":" '{print $1}')
			echo "doing ${num} (${p}) of ${tot}"
			grep -w "${p}" ${blastout}_with_tax | sort -k8 -nr > temp_seq #get all the identities for the hits for a single sequence
			n=$(cut -f7 temp_seq | uniq -c | head -n1 | awk '{print $1}') ##how many of the top identity% are there? If one, pull that and call the hit. If more than one, then go to next.
			top_score=$(cut -f7 temp_seq | head -n1 | awk '{print $1}')
			awk -v d=$top_score '( $8 = d )'  temp_seq | awk -v OFS='\t' '{print $1,$2,$3,$4,$5" "$6,$7,$8}'> temp_choose2
			top_id2=$( sort -k3 -nr temp_choose2 | cut -f3 | head -n1 | awk '{print $1}')
			echo "there are $n hit(s) at best score for sequence ${p}. The top identity % at highest score is $top_id2"
			top_id=$( echo $top_id2 | awk -F"." '{print $1}')
			if [[ $n -eq 1 ]] && [[ $top_id -ge $cutoff ]]
			then	
				ln="head -n1 temp_seq"
				lin=$($ln)
				st=$(echo $lin | awk -v OFS='\t' '{print $1,$3,$5" "$6,$7}')
				taxid=$(echo $lin | awk '{print $7}')
				tax=$(grep -w "^${taxid}" ${ncbi}_r | awk -v OFS='\t' '{print $2,$3,$4,$5,$6}')
				echo "${st}	${tax}	${top_score}	${top_score}" >> ${prefix}_${gene}_best_blast_hits.out
				rm temp_seq
			elif [[ $n -gt 1 ]] && [[ $top_id -ge $cutoff ]]
			then
				echo "more than one equally-good blast hit above ${cutoff}% identity available. Walking up the taxonomy tree to reach consensus."
				id=$(cut -f3 temp_choose2 | sort -nr |  head -n1 | awk '{print $1}')
				grep "$id" temp_choose2 > temp_choose ##pull all the matches with the top identity # for a sequence into a new file##
				spec_number=$(cut -f5 temp_choose | sort | uniq | wc -l | awk '{print $1}') ##how many species are in the equally good hits? If one, pull that and call the hit. if more than one, then go to next.##
				if [[ $spec_number -eq 1 ]]
				then
					spec=$(head -n1 temp_choose |cut -f5 | awk '{print $1,$2}') 
					echo "one species amongst best hits. species is ${spec}."
					ln="head -n1 temp_choose"
					lin=$($ln)
					st=$(echo $lin | awk -v OFS='\t' '{print $1,$3,$5" "$6,$7}')
					taxid=$(echo $lin | awk '{print $7}')
					tax=$(grep -w "^${taxid}" ${ncbi}_r | awk -v OFS='\t' '{print $2,$3,$4,$5,$6}')	
					echo "${st}	${tax}	${top_score}	${top_score}" >> ${prefix}_${gene}_best_blast_hits.out
				else #if multiple species are present in best hits, go through each one and add the taxonomy. This is slow, and is like this because adding a taxmap to the ref database simply does not seem to work.#
					while read z;
					do
						taxid=$(echo ${z} | awk '{print $7}')
						grep -w "^${taxid}" ${ncbi}_r | cut -f2-6 >> temp_ids
					done < temp_choose
					paste temp_choose temp_ids > temp_tax
					rm temp_choose temp_ids
					gen_number=$( cut -f12 temp_tax | sort | uniq | wc -l)
					fam_number=$( cut -f11 temp_tax | sort | uniq | wc -l)
					order_number=$( cut -f10 temp_tax | sort | uniq | wc -l)
					class_number=$( cut -f9 temp_tax | sort | uniq | wc -l)
					phylum_number=$( cut -f8 temp_tax | sort | uniq | wc -l)
				fi
				if [[ ${spec_number} -gt 1 ]] && [[ $gen_number -eq 1 ]]
				then
					gen=$( cut -f12 temp_tax | sort | uniq | awk '{print $1}')
					echo "one genus and ${spec_number} species amongst best hits. taxa is ${gen} sp."
					ln="grep -m1 "$gen" temp_tax"
					lin=$($ln)
					st=$(echo $lin | awk -v OFS='\t' '{print $1,$3}')
					tax=$(echo $lin | awk -v OFS='\t' '{print $9,$10,$11,$12,$13}')
					echo "${st}	${gen} sp.	"NA"	${tax}	${top_score}	${top_score}" >> ${prefix}_${gene}_best_blast_hits.out
				elif [[ ${spec_number} -gt 1 ]] && [[ $fam_number -eq 1 ]]
				then
					fam=$( cut -f11 temp_tax | sort | uniq | awk '{print $1}')
					echo "one family, ${gen_number} genera and ${spec_number} species amongst best hits. taxa is ${fam} sp."
					ln="grep -m1 "$fam" temp_tax"
					lin=$($ln)
					st=$(echo $lin | awk -v OFS='\t' '{print $1,$3}')
					tax=$(echo $lin | awk -v OFS='\t' '{print $9,$10,$11,$12,"NA"}')
					echo "${st}	${fam} sp.	"NA"	${tax}	${top_score}	${top_score}" >> ${prefix}_${gene}_best_blast_hits.out
				elif [[ ${spec_number} -gt 1 ]] && [[ $order_number -eq 1 ]]
				then
					ord=$( cut -f10 temp_tax | sort | uniq | awk '{print $1}')
					echo "one order, ${fam_number} families, ${gen_number} genera and ${spec_number} species amongst best hits. taxa is ${ord} sp."
					ln="grep -m1 "$ord" temp_tax"
					lin=$($ln)
					st=$(echo $lin | awk -v OFS='\t' '{print $1,$3}')
					tax=$(echo $lin | awk -v OFS='\t' '{print $9,$10,$11,"NA","NA"}')
					echo "${st}	${ord} sp.	"NA"	${tax}	${top_score}	${top_score}" >> ${prefix}_${gene}_best_blast_hits.out
				elif [[ ${spec_number} -gt 1 ]] && [[ $class_number -eq 1 ]]
				then
					class=$( cut -f9 temp_tax | sort | uniq | awk '{print $1}')
					echo "one class, ${order_number} orders, ${fam_number} families, ${gen_number} genera and ${spec_number} species amongst best hits. taxa is ${class} sp."
					ln="grep -m1 "$class" temp_tax"
					lin=$($ln)
					st=$(echo $lin | awk -v OFS='\t' '{print $1,$3}')
					tax=$(echo $lin | awk -v OFS='\t' '{print $9,$10,"NA","NA","NA"}')
					echo "${st}	${class} sp.	"NA"	${tax}	${top_score}	${top_score}" >> ${prefix}_${gene}_best_blast_hits.out
				elif [[ ${spec_number} -gt 1 ]] && [[ $phylum_number -eq 1 ]]
				then
					phylum=$( cut -f8 temp_tax | sort | uniq | awk '{print $1}')
					echo "one phylum, ${class_number} classes, ${order_number} orders, ${fam_number} families, ${gen_number} genera and ${spec_number} species amongst best hits. taxa is ${phylum} sp."
					ln="grep -m1 "$phylum" temp_tax"
					lin=$($ln)
					st=$(echo $lin | awk -v OFS='\t' '{print $1,$3}')
					tax=$(echo $lin | awk -v OFS='\t' '{print $9,"NA","NA","NA","NA"}')
					echo "${st}	${phylum} sp.	"NA"	${tax}	${top_score}	${top_score}" >> ${prefix}_${gene}_best_blast_hits.out
				elif [[ ${spec_number} -gt 1 ]] && [[ $phylum_number -gt 1 ]]
				then
					echo "multiple phyla present in equally-good BLAST hits. Designating as no-hit, but reporting best score and top identity."
					st=$( head -n1 temp_tax | awk -v OFS='\t' '{print $1,$3,"No Hit","NA","NA","NA","NA","NA"}')
					echo "${st}	${top_score}	${top_score}" >> ${prefix}_${gene}_best_blast_hits.out
				fi
			elif [[ $top_id -lt $cutoff ]]
			then
				echo "top hit by score has max id% <= ${cutoff}%. Looking for higher %identity hits, potentially at lower scores."
				awk -v d=$cutoff '( $3 > d )' temp_seq > temp_high_pident
				num_hits=$( wc -l temp_high_pident |  awk '{print $1}')
				if [[ $num_hits -eq 0 ]]	
				then
					echo "no BLAST hits of any score above ${cutoff}% identity found."
					if [[ $return_low == "TRUE" ]]
					then
						echo "because you set return_low to TRUE, finding highest score hit regardless of identity."
						id=$(cut -f7 temp_seq | sort -nr |  head -n1 | awk '{print $1}')
						grep -w "${id}$" temp_seq > temp_choose ##pull all the matches with the top score # for a sequence into a new file##
						spec_number=$(cut -f5 temp_choose | sort | uniq | wc -l | awk '{print $1}') ##how many species are in the equally good hits? If one, pull that and call the hit. if more than one, then go to next.##
						if [[ $spec_number -eq 1 ]]
						then
							spec=$(head -n1 temp_choose |cut -f5 | awk '{print $1,$2}') 
							echo "one species amongst best hits. species is ${spec}."
							ln="head -n1 temp_choose"
							lin=$($ln)
							st=$(echo $lin | awk -v OFS='\t' '{print $1,$3,$5" "$6,$7}')
							real_score=$(echo $lin | awk -v OFS='\t' '{print $8}')
							taxid=$(echo $lin | awk '{print $7}')
							tax=$(grep -w "^${taxid}" ${ncbi}_r | awk -v OFS='\t' '{print $2,$3,$4,$5,$6}')	
							echo "${st}	${tax}	${top_score}	${real_score}" >> ${prefix}_${gene}_best_blast_hits.out
						else #if multiple species are present in best hits, go through each one and add the taxonomy. This is slow, and is like this because adding a taxmap to the ref database simply does not seem to work.#
							while read z;
							do
								taxid=$(echo ${z} | awk '{print $7}')
								grep -w "^${taxid}" ${ncbi}_r | cut -f2-6 >> temp_ids
							done < temp_choose
							paste temp_choose temp_ids > temp_tax
							rm temp_choose temp_ids
							gen_number=$( cut -f12 temp_tax | sort | uniq | wc -l)
							fam_number=$( cut -f11 temp_tax | sort | uniq | wc -l)
							order_number=$( cut -f10 temp_tax | sort | uniq | wc -l)
							class_number=$( cut -f9 temp_tax | sort | uniq | wc -l)
							phylum_number=$( cut -f8 temp_tax | sort | uniq | wc -l)
						fi
						if [[ ${spec_number} -gt 1 ]] && [[ $gen_number -eq 1 ]]
						then
							gen=$( cut -f12 temp_tax | sort | uniq | awk '{print $1}')
							echo "one genus and ${spec_number} species amongst best hits. taxa is ${gen} sp."
							ln="grep -m1 "$gen" temp_tax"
							lin=$($ln)
							st=$(echo $lin | awk -v OFS='\t' '{print $1,$3}')
							real_score=$(echo $lin | awk -v OFS='\t' '{print $8}')
							tax=$(echo $lin | awk -v OFS='\t' '{print $9,$10,$11,$12,$13}')
							echo "${st}	${gen} sp.	"NA"	${tax}	${top_score}	${real_score}" >> ${prefix}_${gene}_best_blast_hits.out
						elif [[ ${spec_number} -gt 1 ]] && [[ $fam_number -eq 1 ]]
						then
							fam=$( cut -f11 temp_tax | sort | uniq | awk '{print $1}')
							echo "one family, ${gen_number} genera and ${spec_number} species amongst best hits. taxa is ${fam} sp."
							ln="grep -m1 "$fam" temp_tax"
							lin=$($ln)
							st=$(echo $lin | awk -v OFS='\t' '{print $1,$3}')
							real_score=$(echo $lin | awk -v OFS='\t' '{print $8}')
							tax=$(echo $lin | awk -v OFS='\t' '{print $9,$10,$11,$12,"NA"}')
							echo "${st}	${fam} sp.	"NA"	${tax}	${top_score}	${real_score}" >> ${prefix}_${gene}_best_blast_hits.out
						elif [[ ${spec_number} -gt 1 ]] && [[ $order_number -eq 1 ]]
						then
							ord=$( cut -f10 temp_tax | sort | uniq | awk '{print $1}')
							echo "one order, ${fam_number} families, ${gen_number} genera and ${spec_number} species amongst best hits. taxa is ${ord} sp."
							ln="grep -m1 "$ord" temp_tax"
							lin=$($ln)
							st=$(echo $lin | awk -v OFS='\t' '{print $1,$3}')
							real_score=$(echo $lin | awk -v OFS='\t' '{print $8}')
							tax=$(echo $lin | awk -v OFS='\t' '{print $9,$10,$11,"NA","NA"}')
							echo "${st}	${ord} sp.	"NA"	${tax}	${top_score}	${real_score}" >> ${prefix}_${gene}_best_blast_hits.out
						elif [[ ${spec_number} -gt 1 ]] && [[ $class_number -eq 1 ]]
						then
							class=$( cut -f9 temp_tax | sort | uniq | awk '{print $1}')
							echo "one class, ${order_number} orders, ${fam_number} families, ${gen_number} genera and ${spec_number} species amongst best hits. taxa is ${class} sp."
							ln="grep -m1 "$class" temp_tax"
							lin=$($ln)
							st=$(echo $lin | awk -v OFS='\t' '{print $1,$3}')
							tax=$(echo $lin | awk -v OFS='\t' '{print $9,$10,"NA","NA","NA"}')
							real_score=$(echo $lin | awk -v OFS='\t' '{print $8}')
							echo "${st}	${class} sp.	"NA"	${tax}	${top_score}	${ereal_score}" >> ${prefix}_${gene}_best_blast_hits.out
						elif [[ ${spec_number} -gt 1 ]] && [[ $phylum_number -eq 1 ]]
						then
							phylum=$( cut -f8 temp_tax | sort | uniq | awk '{print $1}')
							echo "one phylum, ${class_number} classes, ${order_number} orders, ${fam_number} families, ${gen_number} genera and ${spec_number} species amongst best hits. taxa is ${phylum} sp."
							ln="grep -m1 "$phylum" temp_tax"
							lin=$($ln)
							st=$(echo $lin | awk -v OFS='\t' '{print $1,$3}')
							tax=$(echo $lin | awk -v OFS='\t' '{print $9,"NA","NA","NA","NA"}')
							real_score=$(echo $lin | awk -v OFS='\t' '{print $8}')
							echo "${st}	${phylum} sp.	"NA"	${tax}	${top_score}	${real_score}" >> ${prefix}_${gene}_best_blast_hits.out
						elif [[ ${spec_number} -gt 1 ]] && [[ $phylum_number -gt 1 ]]
						then
							echo "multiple phyla present in equally-good BLAST hits. Designating as no-hit, but reporting best score and top identity."
							st=$( head -n1 temp_tax | awk -v OFS='\t' '{print $1,$3,"No Hit","NA","NA","NA","NA","NA"}')
							echo "${st}	${top_score}	${top_score}" >> ${prefix}_${gene}_best_blast_hits.out
							real_score=$(echo $lin | awk -v OFS='\t' '{print $8}')
						fi
					else
						echo "because you set return_low to FALSE, designating no-hit, but reporting best score and highest pident."
						st=$( head -n1 temp_seq | awk -v OFS='\t' '{print $1,$3,"No Hit","NA","NA","NA","NA","NA","NA"}')
						echo "${st}	${top_score}	0" >> ${prefix}_${gene}_best_blast_hits.out	
					fi
				else
					echo "found BLAST hit(s) with ${cutoff}% id or higher, but at below highest score. Reporting both real and top score along with identity."
					id=$(cut -f3 temp_seq | sort -nr | head -n1 | awk '{print $1}')
					grep "${id}	" temp_seq > temp_choose ##pull all the matches with the top score # for a sequence into a new file##
					spec_number=$(cut -f5 temp_choose | sort | uniq | wc -l | awk '{print $1}') ##how many species are in the equally good hits? If one, pull that and call the hit. if more than one, then go to next.##
					if [[ $spec_number -eq 1 ]]
					then
						spec=$(head -n1 temp_choose |cut -f5 | awk '{print $1,$2}') 
						echo "one species amongst best hits. species is ${spec}."
						ln="head -n1 temp_choose"
						lin=$($ln)
						st=$(echo $lin | awk -v OFS='\t' '{print $1,$3,$5" "$6,$7}')
						real_score=$(echo $lin | awk -v OFS='\t' '{print $8}')
						taxid=$(echo $lin | awk '{print $7}')
						tax=$(grep -w "^${taxid}" ${ncbi}_r | awk -v OFS='\t' '{print $2,$3,$4,$5,$6}')	
						echo "${st}	${tax}	${top_score}	${real_score}" >> ${prefix}_${gene}_best_blast_hits.out
					else #if multiple species are present in best hits, go through each one and add the taxonomy. This is slow, and is like this because adding a taxmap to the ref database simply does not seem to work.#
						while read z;
						do
							taxid=$(echo ${z} | awk '{print $7}')
							grep -w "^${taxid}" ${ncbi}_r | cut -f2-6 >> temp_ids
						done < temp_choose
						paste temp_choose temp_ids > temp_tax
						rm temp_choose temp_ids
						gen_number=$( cut -f12 temp_tax | sort | uniq | wc -l)
						fam_number=$( cut -f11 temp_tax | sort | uniq | wc -l)
						order_number=$( cut -f10 temp_tax | sort | uniq | wc -l)
						class_number=$( cut -f9 temp_tax | sort | uniq | wc -l)
						phylum_number=$( cut -f8 temp_tax | sort | uniq | wc -l)
					fi
					if [[ ${spec_number} -gt 1 ]] && [[ $gen_number -eq 1 ]]
					then
						gen=$( cut -f12 temp_tax | sort | uniq | awk '{print $1}')
						echo "one genus and ${spec_number} species amongst best hits. taxa is ${gen} sp."
						ln="grep -m1 "$gen" temp_tax"
						lin=$($ln)
						st=$(echo $lin | awk -v OFS='\t' '{print $1,$3}')
						real_score=$(echo $lin | awk -v OFS='\t' '{print $8}')
						tax=$(echo $lin | awk -v OFS='\t' '{print $9,$10,$11,$12,$13}')
						echo "${st}	${gen} sp.	"NA"	${tax}	${top_score}	${real_score}" >> ${prefix}_${gene}_best_blast_hits.out
					elif [[ ${spec_number} -gt 1 ]] && [[ $fam_number -eq 1 ]]
					then
						fam=$( cut -f11 temp_tax | sort | uniq | awk '{print $1}')	
						echo "one family, ${gen_number} genera and ${spec_number} species amongst best hits. taxa is ${fam} sp."
						ln="grep -m1 "$fam" temp_tax"
						lin=$($ln)
						st=$(echo $lin | awk -v OFS='\t' '{print $1,$3}')
						real_score=$(echo $lin | awk -v OFS='\t' '{print $8}')
						tax=$(echo $lin | awk -v OFS='\t' '{print $9,$10,$11,$12,"NA"}')
						echo "${st}	${fam} sp.	"NA"	${tax}	${top_score}	${real_score}" >> ${prefix}_${gene}_best_blast_hits.out
					elif [[ ${spec_number} -gt 1 ]] && [[ $order_number -eq 1 ]]
					then
						ord=$( cut -f10 temp_tax | sort | uniq | awk '{print $1}')
						echo "one order, ${fam_number} families, ${gen_number} genera and ${spec_number} species amongst best hits. taxa is ${ord} sp."
						ln="grep -m1 "$ord" temp_tax"
						lin=$($ln)
						st=$(echo $lin | awk -v OFS='\t' '{print $1,$3}')
						real_score=$(echo $lin | awk -v OFS='\t' '{print $8}')
						tax=$(echo $lin | awk -v OFS='\t' '{print $9,$10,$11,"NA","NA"}')
						echo "${st}	${ord} sp.	"NA"	${tax}	${top_score}	${real_score}" >> ${prefix}_${gene}_best_blast_hits.out
					elif [[ ${spec_number} -gt 1 ]] && [[ $class_number -eq 1 ]]
					then
						class=$( cut -f9 temp_tax | sort | uniq | awk '{print $1}')
						echo "one class, ${order_number} orders, ${fam_number} families, ${gen_number} genera and ${spec_number} species amongst best hits. taxa is ${class} sp."
						ln="grep -m1 "$class" temp_tax"
						lin=$($ln)
						st=$(echo $lin | awk -v OFS='\t' '{print $1,$3}')
						tax=$(echo $lin | awk -v OFS='\t' '{print $9,$10,"NA","NA","NA"}')
						real_score=$(echo $lin | awk -v OFS='\t' '{print $8}')
						echo "${st}	${class} sp.	"NA"	${tax}	${top_score}	${real_score}" >> ${prefix}_${gene}_best_blast_hits.out
					elif [[ ${spec_number} -gt 1 ]] && [[ $phylum_number -eq 1 ]]
					then
						phylum=$( cut -f8 temp_tax | sort | uniq | awk '{print $1}')
						echo "one phylum, ${class_number} classes, ${order_number} orders, ${fam_number} families, ${gen_number} genera and ${spec_number} species amongst best hits. taxa is ${phylum} sp."
						ln="grep -m1 "$phylum" temp_tax"
						lin=$($ln)
						st=$(echo $lin | awk -v OFS='\t' '{print $1,$3}')
						tax=$(echo $lin | awk -v OFS='\t' '{print $9,"NA","NA","NA","NA"}')
						real_score=$(echo $lin | awk -v OFS='\t' '{print $8}')
						echo "${st}	${phylum} sp.	"NA"	${tax}	${top_score}	${real_score}" >> ${prefix}_${gene}_best_blast_hits.out
					elif [[ ${spec_number} -gt 1 ]] && [[ $phylum_number -gt 1 ]]
					then
						echo "multiple phyla present in equally-good BLAST hits. Designating as no-hit, but reporting best score and top identity."
						st=$( head -n1 temp_tax | awk -v OFS='\t' '{print $1,$3,"No Hit","NA","NA","NA","NA","NA"}')
						echo "${st}	${top_score}	${top_score}" >> ${prefix}_${gene}_best_blast_hits.out
						real_score=$(echo $lin | awk -v OFS='\t' '{print $8}')
					fi
				fi	

			fi
		fi
	done < temp_seqlist
	rm temp_seqlist
	
 	echo "done with choosing best blast hits, now creating and formatting outfiles."
	h1="head -n1 ${prefix}_${gene}_best_blast_hits.out"
	head1=$($h1)
	tail -n +2 ${prefix}_${gene}_best_blast_hits.out > ${prefix}_${gene}_best_blast_hits.out.r
	echo "sequence	${head1}" > ${prefix}_${gene}_best_blast_hits.txt
	while read c;
	do
		seqnum=$( echo ${c} | awk '{print $1}') 
		seq=$( grep -w -A1 ">${seqnum}" ${prefix}_${gene}_combined_ASVs.fasta | tail -n1 | awk '{print $1}')
		echo "${seq}	${c}" >> ${prefix}_${gene}_best_blast_hits.txt
	done < ${prefix}_${gene}_best_blast_hits.out.r
	rm temp_sq ${prefix}_${gene}_best_blast_hits.out ${prefix}_${gene}_best_blast_hits.out.r
	
 	echo "sample	sequence	reads	identity	taxa	taxid	top_score	chosen_score	phylum	class	order	family	genus" >  ${prefix}_${gene}_taxatable.txt
	for fil in *_F_filt.fastq_seqs.txt
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
			taxline=$(grep -w  -m1 "$seq" ${prefix}_${gene}_best_blast_hits.txt | awk -v OFS='\t' '{print $3,$4" "$5,$6,$7,$8,$9,$10,$11,$12,$13}')
			echo "$base	$seq	$reads	$taxline" >> ${prefix}_${gene}_taxatable.txt
			x=$(( $x + 1 ))
		done
	done
	echo "done with ${gene}. You can find taxa tables and raw ASV tables in your project directory/reports."
	rm temp*
	mkdir ${dirr}/${gene}_dada_out/sample_seqfiles
	mv *_seqs.txt ${dirr}/${gene}_dada_out/sample_seqfiles
	cd ${dirr}
	cp ${dirr}/${gene}_dada_out/*_taxatable.txt ${dirr}/results_tables/

done < list_of_genes.txt
