#!/bin/bash

###This is the jankiest script of the whole bunch###
#BEFORE RUNNING THIS SCRIPT CHECK YOUR RESULTS FILE. COUNT HOW MANY LINES BETWEEN THE "Query=" line and the top hit!!!!!#
#If that number is NOT 6, for instance, if it is 8 (it depends on local/remote and sometimes just how NCBI feels that day), you need to change the "line=$(( $lnno + 6 ))" below from 6 to your # of lines.
#This script needs heavy modifications before it can be used blithely.



fasta=your_project_seqs_only_formatted.fasta ##this is the fasta you input into BLAST. It has your query names in it.
results=your_project_blast_results.out ##this is your blast results file
out=your_project_with_taxa.txt ##this is the output of this file: a list of sequences and their associated taxa, and the % identity of the hit.
echo "sequence  taxa    identity" > ${out} ##this is the format of the output. Three tab-delimited columns with sequence, taxa, and identity.
  y=1
        seqy=2
        recy=2
        z=$(wc -l $fasta | awk '{print $1}')
        while [[ $y -le $z ]] ##this is a loop that reads through the fasta file and for each sequence, finds the results and assesses them##
        do
		sqlin="sed -n ${y}p $fasta" #grab the yth line of the fasta file. This is a counter (1 to n) that skips every other line, because we only want the headers of the sequences.
                que=$($sqlin) #turn the fasta file line into a variable we can manipulate
                query="${que:1}" #chomp the ">" off the front
                echo "query is ${query}" #echo the query so we know we're grabbing the right thing. This is a debug line.
		lnno=$( grep -n -m1 "${query}" $results | awk -F":" '{print $1}') #we want the line number of the results file that has our query's results in it. I tried a lot of ways of doing this but this is easiest: grep for the query, with the -n flag to grab the line number and -m1 to only grab the first match. Cut the line number (it is separated by a ":" from the match itself). Store that as line number
                echo "line number is $lnno" #debug line, list the linenumber.
                if gsed -n "/Query= ${query}$/, /Query= /p" $results | grep -q "No hit"; #check if there was a hit. For local blast, sometimes you get no hit. If no hit, just say "no hit" and leave it out of the final thing. Although I might modify this so that the results just say no hit, which would simplify some of my downstream stuff.
		then
		echo "no hit for ${query}"
		y=$(( $y + 2 )) #increase the counters for the loop
                recy=$(( $recy + 1 ))
                seqy=$(( $seqy + 2))
	else
		line=$(( $lnno + 6 )) ##this is the line to change and check!! this is important. In general, this value should be 6 for local blast and 8 for remote blast. This finds the top blast hit for each query by moving n lines down from the "Query:" line of the results file. This is the janky bit but it's very difficult to figure out another way.
                ans="sed -n ${line}p ${results}" #yank the best hit line#
                res=$($ans)
                taxa=$(echo $res | awk -F" " '{print $1}' | awk -F"_" '{print $2,$3}') #this separates the best hit line by spaces and pulls the first two words as the taxname. This is a good, but not perfect, method of doing it.
		echo "$taxa" #debug line, list the taxa
                sq="sed -n ${seqy}p $fasta" #this grabs the sequence from the fasta file
                seq=$($sq)
                #ct="sed -n ${recy}p ${records}"
		identity=$(gsed -n "/Query= ${query}$/, /Query= /p" $results | grep -m1 "Identities"  | awk -F" " '{print $4}' | cut -c2- | awk -F")" '{print $1}') #This grabs the associated identity from the top hit. This is also janky, because it only works if you ask BLAST to return at least one alignment (you can just ask BLAST to return descriptions instead). But it does at least usually, mostly work. However, depending on remote/local, sometimes this line needs to be changed to include or remove a space after the "Query=${query}". So, also janky.
                echo "$identity"
                echo "${seq}       ${taxa}      ${identity}" >> ${out} ##print to your output. Sequence, taxa, identity.
                y=$(( $y + 2 )) ##increase all your counters##
                recy=$(( $recy + 1 ))
                seqy=$(( $seqy + 2))
		fi
	done
