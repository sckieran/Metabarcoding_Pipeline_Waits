# Metabarcoding_Pipeline_Waits
A rough pipeline and tutorial to analyze metabarcoding data using dada2. This explains my system and can help someone parse the scripts but is not a true, workable pipeline.

This pipeline is designed to run in the terminal on Mac OS Monterey. It will need to be considerably modified for each project.

Software needs: fastq-multx, bbmap, gsed, command-line BLAST. R packages: tidyverse, dada2, rentrez,tidyr, stringr,ggarrange

**Step 0: Create a reference database from sequences available on NCBI.**
        **0.1 Script: query_rentrez.R** and its associated project-specific files.
            **Input:** List of taxa (taxids optional but preferred), list of genes with "gene terms". This input looks like a list of columns. Each column corresponds to a gene, each row is different ways to describe that gene. This just provides the best search radius. So 12S would include "12S" "12S Ribosomal RNA" "12S Mitochondrial Sequence" etc etc. COI is the worst for this.
         ** Output:**  Up to 10 fasta files for every taxa on the taxalist with hits in NCBI for the gene in question, plus summary files of which taxa had sequences, how many, etc. Can be easily modified to just return summaries or just return fastas. Runs faster if you add an API key. Sometimes pings NCBI too fast without it. For broadest use, will probably use API key for personal script but leave it out for general script and just add a sleep 0.1 command to keep it from hitting the rate limit of 3/sec. Taxa fastas can fill up a folder quickly so it's a good idea to make one specially for that and direct your script there.
       ** 0.2** cat taxa files together. (cat \*.fasta > all_refseqs.fasta)
       ** 0.3** build NCBI local databse.
        Syntax is makeblastdb -in all_refseqs.fasta -dbtype nucl -parse_seqids -out your_reference_name
        if you use the query_rentrez script for this, the sequences are already perfectly formatted for this.
        
   **Step 1: Demultiplex.** No script, uses demux-by-name.sh (bbmap) or fastq-multx. I write this one separately for each sequencing run, since it's one line of code specific to the project.

**Step 2: Rename, if bbmap**.
        2.1: Script: rename.sh. Input: List of barcodes and associated samples. Output: Files are renamed from barcode to sample name. Only need this step if demult with demux-by-name instead of fastq-multx. Also need to use this step to make sure all files have unique names - if multiple plates, all PCRN/P have to have unique identifiers. Also need to move anything you don't want to analyze out of your folder full of sequences - I have an "unused" folder where I keep the unmatched.fastqs from demultiplexing, but they could also be removed.

**Step 3: DADA processing and BLAST prep. **
  ** 3.1 Script:** DADA2_new.R.** I make one of these for each project, because there is too much to change to do it on the command line on the fly. Also, run this interactively. Lots of intermediate graphs.
                **Input:** folder with uniquely-named F/R files, and the pattern of their names (for me, generally filename_R1.fastq/filename_R2.fastq).
                **Output**: A comma-separated list of unique ASVs passing filter for each sample (merged F/Rs), along with their read counts. A list of all samples passing a looser filter, one per line. A tab-separated table of read counts retained and lost at each step. A list of sequences that were filtered out based on the read count filter.
                Notes: This is the big processing script. It simply isn't written to be run in a big lump. Even in the DADA tutorial, you're supposed to look at the output of these things, look at the tracking table, and consider what your data needs in terms of filters. Locally, this script will top out my Macbook Pro if I have more than 6 or 7 plates (600-700 samples). If we want this to be a spoon-fed pipeline, we need to think about cutting it up into smaller chunks with the filtering outputs delineating the ends of each script. Also, this script needs to be "variable-ized", because right now large chunks of the output names need to be changed by hand throughout the script. So we need to set some global variables at the start.
       ** 3.2 Script: reformat_for_BLAST.sh** This makes a list of all the unique sequences present in at least one sample and formats them as a fasta file for blasting.
                Input: The \*_seqs.txt output of the DADA2 script. Only works if they're named like that, or you can change the name pattern. The script will grab any file that matches the name pattern, so be careful about what else you put in there.
                Output: A fasta-formatted file of sequences to blast. The sequences will be named numerically, so they will be >seq_1, >seq_2....>seq_n.

**Step 4: BLAST**
       ** 4.1 Script: blast.sh.** I wrote one of these but I don't actually usually use them. It's a single line and it depends on whether you're doing local or remote blast. syntax is blastn -db [database] -query [fasta] -out [results file] -num_alignments=5 -num_descriptions=5 [-remote]. You can change the num_alignments but I recommend >= 1. You can change the number of descriptions, 1 will pull only the top hit. BLAST warns you that it should be >1. Lots of filters you can add but I like to be broad and filter later. Output is a horribly-formatted human-readable results file that has a large header and then lists each query, the top hits, and their alignments.

**Step 5: **
Make a taxa table. This script needs significant modificants to manage edge cases. If you do local blast, line  A lot of this step is really janky due to how weird the BLAST output results file is and how difficult it is to parse. I am thinking about a total re-work of these scripts but I'm in the very early stages.
        5.1 Script: add_taxa_to_blast_results.sh
                Inputs: The fasta that was the input to the BLAST search and the BLAST results file.
                Output: A tab-separated file with a list of unique sequences, the taxa they blast to (top blast hit), and the identity % of the match. This is the jankiest script because the BLAST results file changes slightly based on the format of the fasta input and whether it's remote or local, and a few other parameters, like the weather, and how BLAST feels that day, and whether you use vim or emacs etc etc. The script is extensively commented but still terrible.
        5.2 make_taxa_table.sh
                Inputs: _seqs.txt files for each sample (the output of the DADA2 script), and the blast-results-with-taxonomy file that is the output of step 5.1
                Outputs: a tab-separated table with columns for replicate, sequence, reads, taxa and identity. Every sequence from every replicate that passed the DADA filters is included in this table, with the read counts. Still finding the best way to handle no-hits but I've got some ideas and it works ok for now.

Step 6: Analyses
        Step 6.1 I usually open that taxa table in excel and add common names (HR taxa, "human-readable") and other data by hand. I literally do common names by making a unique list of taxa in a new worksheet and just googling. Excel has a nice in-house googling function that I use sometimes (CMD-CTRL-L), but it can't fill cells for you (it just returns the top google hit for the contest of the cell in a sidebar for you to look at). For some stuff I also add additional taxa info, like family (for plants), or prey type (for neil's stuff, bird v mammal v PCRP). If I'm feeling fancy, I add a column for colors for R analysis later. I also add a "sample" column, which is just replicate-agnostic. I usually do some pivot tables with read counts and taxa. Once I feel good about it, I send it off, and usually save a tab-delimited text version for loading into R.


