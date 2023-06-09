#!/bin/bash

##this is remote blast. For local blast, change -db to the dbname you set when you made your blastdb. If you didn't assign a name, it's the name of the fasta, including the '.fasta' part.##


##queries blast for each sequence. Must have ncbi command line blast installed. Pulls the top 5 hits and shows the alignment of the best one.
blastn -db nt -query your_project_seqs_only_formatted.fasta -out your_project_blast_results.out -num_alignments=1 -num_descriptions=5 -remote
