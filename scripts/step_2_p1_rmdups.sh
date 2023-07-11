#!/bin/bash


prefix=$1
gene=$2

filename=${prefix}_${gene}_sequence_database.fasta

grep ">" $filename | cut -d" " -f 1-5 | sort | uniq -c | awk '($1 > 1)' | awk -F" " '{print $2}' > dup_fs

nl=$( wc -l dup_fs | awk '{print $1}' )
if [[ $nl -eq 0 ]]
then
echo "checked for duplicate fasta entries for ${gene}. No duplicates found. Making blast databases."
else
dfile=dup_fs
grep -n ">" $filename > n_lines
line_file=n_lines

while read p;
do
acc=${p}
grno=$(grep -c "$acc" $line_file | awk '{print $1}')
y=1
while [[ $y -lt $grno ]]
do
        start_line=$(grep -m${y} "$acc" $line_file | awk -F":" '{print $1}' | tail -n1)
        end1=$(grep -m${y} -A1 "$acc" $line_file | grep -v "$acc" |  awk -F":" '{print $1}' | tail -n1)
        end_line=$(( $end1 - 1 ))
        echo "start line is $start_line and end line is $end_line"
        sed -e "${start_line},${end_line}d" $filename > temp_file
        grep -n ">" temp_file > $line_file
        mv temp_file $filename
        y=$(( $y + 1 ))
done

done < $dfile
fi
