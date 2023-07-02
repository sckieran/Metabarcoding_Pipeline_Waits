#!/bin/bash
meta=$1

project_dir=$(sed -n 1p $meta | awk -F"=" '{print $2}' | awk -F" #" '{print $1}')
rra_cutoff=$(sed -n 2p $meta | awk -F"=" '{print $2}' | awk -F" #" '{print $1}')  #relative read abundance cutoff. Default is 0.001 (0.1%). That means if a SAMPLE has 10,000 reads total, this removes any ASV with <10 reads. It is reasonably strict if you have high read counts.
filt_dir==$(sed -n 3p $meta | awk -F"=" '{print $2}' | awk -F" #" '{print $1}') #path to filtered fastas will be /path/to/your/dir/your_gene/filtered/, path to outfiles will be /path/to/your/dir/your_gene_dada_out/
truncLenF=$(sed -n 4p $meta | awk -F"=" '{print $2}' | awk -F" #" '{print $1}') #truncation length for R1s. change this based on your run length (PE 75, 150 or 300), your overlap (read length - amplicon length) and your quality report output
truncLenR=$(sed -n 5p $meta | awk -F"=" '{print $2}' | awk -F" #" '{print $1}') #truncation length for R2s. see above
maxEEF=$(sed -n 6p $meta | awk -F"=" '{print $2}' | awk -F" #" '{print $1}') #expected errors for R1s. Default is 2. See DADA2 tutorial for more detailed filtering instructions: https://benjjneb.github.io/dada2/tutorial_1_8.html
maxEER=$(sed -n 7p $meta | awk -F"=" '{print $2}' | awk -F" #" '{print $1}') #expected errors for R2s
truncQ=$(sed -n 8p $meta | awk -F"=" '{print $2}' | awk -F" #" '{print $1}') #truncate read at first base with quality score Q
rmphix=$(sed -n 9p $meta | awk -F"=" '{print $2}' | awk -F" #" '{print $1}') #remove reads that match the phix genome, set to TRUE if you had your sequence spike in PhiX, common for amplicon sequencing to increase complexity
compress=$(sed -n 10p $meta | awk -F"=" '{print $2}' | awk -F" #" '{print $1}') #set to TRUE if you want to compress outfiles
multithread=$(sed -n 11p $meta | awk -F"=" '{print $2}' | awk -F" #" '{print $1}') #set to TRUE if you want to multithread, generally true on clusters and macs, false on windows.
trimLeft=$(sed -n 12p $meta | awk -F"=" '{print $2}' | awk -F" #" '{print $1}') #From dada: Default 0. The number of nucleotides to remove from the start of each read. If both truncLen and trimLeft are provided, filtered reads will have length truncLen-trimLeft.
trimRight=$(sed -n 13p $meta | awk -F"=" '{print $2}' | awk -F" #" '{print $1}') #From dada: Default 0. The number of nucleotides to remove from the end of each read. If both truncLen and trimRight are provided, truncation will be performed after trimRight is enforced.
maxLen=$(sed -n 14p $meta | awk -F"=" '{print $2}' | awk -F" #" '{print $1}') #From dada: Default Inf (no maximum). Remove reads with length greater than maxLen. maxLen is enforced before trimming and truncation.
minLen=$(sed -n 15p $meta | awk -F"=" '{print $2}' | awk -F" #" '{print $1}') #(Optional). Default 20. Remove reads with length less than minLen. minLen is enforced after trimming and truncation.


touch filtercomm
echo "out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN=0," > filtercomm


if [[ -z ${truncLenF} ]] && [[ -z ${truncLenR} ]]
then
	echo "both trunlenF and trunlenR are empty. Omitting trunlen from filtering parameters."
elif [[ -z ${truncLenF} ]] && [[ ! -z ${truncLenR} ]]
then	
	echo "only trunlenR supplied. Not truncating R1s."
	echo "truncLen=c(0,${trunLenR})," >> filtercomm
elif [[ ! -z ${truncLenF} ]] && [[  -z ${truncLenR} ]]
then	
	echo "only trunlenF supplied, not truncating R2s."
	echo "truncLen=c(${trunclenF},0)," >> filtercomm
fi
if [[ -z ${maxEEF} ]] && [[ -z ${maxEER} ]]
then
        echo "both maxEEF and maxEER are empty. Omitting maxEE from filtering parameters."
elif [[ -z ${maxEEF} ]] && [[ ! -z ${maxEER} ]]
then
        echo "maxEE=c(0,${maxEER})," >> filtercomm
elif [[ ! -z ${maxEEF} ]] && [[ -z ${maxEER} ]]
then    
    echo "truncLen=c(${maxEEF},0)," >> filtercomm
fi

if [[ ! -z ${truncQ} ]]
then
echo "truncQ=${truncQ}," >> filtercomm
fi

if [[ ! -z ${rmphix} ]]
then
echo "rm.phix=${rmphix}," >> filtercomm
fi

if [[ ! -z ${compress} ]]
then
echo "compress=${compress}," >> filtercomm
fi

if [[ ! -z ${multithread} ]]
then
echo "multithread=${multithread}," >> filtercomm
fi

if [[ ! -z ${trimLeft} ]]
then
echo "trimLeft=${trimLeft}," >> filtercomm
fi

if [[ ! -z ${trimRight} ]]
then
echo "trimRight=${trimRight}," >> filtercomm
fi

if [[ ! -z ${maxLen} ]]
then
echo "maxLen=${maxLen}," >> filtercomm
fi

if [[ ! -z ${minLen} ]]
then
echo "minLen=${minLen}," >> filtercomm
fi

tr '\n' ' ' < filtercomm > filtline 
sed -i 's/, $/)/' filtline	
