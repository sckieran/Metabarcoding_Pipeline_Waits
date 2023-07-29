# Waits Lab Metabarcoding Training Pipeline
Shannon Blair 2023

## Overview of Pipeline

A pipeline and tutorial to analyze metabarcoding data using pears and fastx-collapser with a limited local reference database. It is intended to be an analysis guide for students working on their own projects. It currently only runs on the UI RCDS cluster.


This pipeline assumes reasonably good resolution of locus data and is designed as a "first pass" at your metabarcoding data, not a finished product. The way this pipeline assesses BLAST hits is both **naive and conservative**. If more than one equally-likely BLAST hit is available (judged by bitscore), it will attempt to assign taxonomy to each hit using the output of ncbitax2lin, then walk up the taxonomy tree until it finds a consensus rank. If no consensus is available at the phylum level, it calls a "no-hit". Therefore **this pipeline only works with reference libraries made using this pipeline.** However, you can use your own FASTAs as input to step 2 of this pipeline if you want to skip querying rentrez, as long as each sequence in that FASTA begins with a header line in the format >[Valid NCBI Accession Number]. The pipeline also returns all the equally-good top hits for each sequence, so you can examine them individually and assess the validity of the classifier.


# FAQ
**What if I want my database to include every organism available for X gene?**

This pipeline is for limited-taxa reference databases or remote querying of NCBI. It is not designed to manage the curation of a large database that includes, say, all inverts with COI sequences in NCBI. If you want a larger database, consider using another reference database building tool, for instance [RESCRIPt](https://github.com/bokulich-lab/RESCRIPt) or [bcDatabaser](https://bcdatabaser.molecular.eco/). This pipeline is being actively developed to better improve performance with other databases. Any FASTA file that contains unique sequences each with a valid NCBI accession will work in step 2 of the pipeline, but runtimes for step 2 will be very long for huge databases (>1000 taxa).


**What if I have a few of my own in-house sequences I want to add to my database?**

That's fine, but it will add some tedioius work upfront. First, check if the species you are including have taxIDS assigned in NCBI. All taxa with at least one sequence in the NCBI nucleotide or sra database have a taxid, and you can look it up here: [here](https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi). Put your homebrewed sequences in a file called extras_gene1_sequences.fasta with the header format ">unique_identifier species name taxid=unique_taxid" and put it in your reference_database directory. If your species has a taxID assigned in NCBI, use that as the unique_taxid. Otherwise, you can assign one yourself, I recommend a series of letters and numbers separated by underscores. Do not use a strictly numeric value, you will almost definitely accidentally assign a real taxID to your species. You can assign any unique_identifier to the sequence (to replace the accession number), but it must not include spaces or special characters and each sequence must be uniquely named. Next, also in your reference database directory, install and run ncbitax2lin. Then, unzip and modify the ncbitax2lin file. Add the unique_taxids you assigned to the end of the file and add taxonomic information for each species following the comma-separated format of the ncbitax2lin file (described in the header of the file). This is tedious to do by hand, so we really recommend submitting your homebrew sequences to NCBI if they're long enough!

**What if I want to use remote BLAST to query all of NCBI?**

The pipeline does not currently do remote blast. Because there is no pre-BLAST filtering in this version of the pipeline, many hundreds of thousands of ASVs may be returned from a large dataset. It is not currently feasible to remote-BLAST more than ~3,000 sequences at once due to the time it takes to ping NCBI's servers, and the likelihood of NCBI throwing an error and cancelling the job.

**Why doesn't this pipeline use e-value instead of bitscore or pident?**
E-values help asses likelihood of homology and are affected by database size. The conserved regions of metabarcoding primers almost guarantee homology above any standard e-value cutoff, and these databases are generally small, so e-values are high across the board.

**I want to do 16S bacterial metagenomics. Is this pipeline right for me?** Probably not. There is an immense, incredibly well-maintained array of resources for 16S microbial amplicon sequencing, some commercial and some open source, that will be infinitely better than this pipeline.


## Inputs and Preparation

### Inputs
This pipeline requires the following inputs:
- Demultiplexed, appropriately-named **paired-end metabarcoding data in fastq (or fastq.gz) format**, with adapters trimmed. Each gene or primer set should be analyzed separately, and should be in a different folder called gene1, gene2...geneN.
- The output of a single run of **ncbitax2lin**. This script will attempt to install and run ncbitax2lin if it can't find a taxonomy file. If you are having issues with permissions and install, you can run it yourself by copying the commands from [this readme](https://github.com/zyxue/ncbitax2lin), the script assumes a name of "ncbi_lineages_[date_of_utcnow].csv.gz" but is agnostic to the date.
- A **taxa file** containing the scientific names of your target taxa, one per line, with the header "taxname", if you want to run step 1 and query NCBI for your fasta.
- A **"genelist" file** with lists of common terms for your target genes. One gene per column, as many permutations on the gene as you'd like, one per line (ie, "Cytochrome Oxidase I", "COI", "COX1"). See the example files for a template. The columns should have a header (which you can repeat in the body) that is a short, human-readable name of the gene that contains no spaces, slashes, quote marks or other special characters. For example, instead of heading your column "Cytochrome Oxidase I", head it "COI". 

Data is often demultiplexed by the sequencing service company at no (or minor) cost. However, if your data has not been demultiplexed, we recommend using either [fastq-multx](https://github.com/brwnj/fastq-multx) or the demux-by-name function of [BBMap](https://github.com/BioInfoTools/BBMap). Some demultiplexers (fastq-multx, for instance) do automatic trimming, others do not, which is why this pipeline only accepts pre-trimmed data. We recommend using [cutadapt](https://cutadapt.readthedocs.io/en/stable/) to trim. 

### Software and Installation

You will need
- Command-line git (pre-installed on UI RCDS Cluster, UCD FARM and Barbera)
- Command-Line BLAST [Installation Instructions Here](https://www.ncbi.nlm.nih.gov/books/NBK279690/) (pre-installed on UI RCDS Cluster, UCD FARM and Barbera)
- NCBItax2lin: [see here](https://github.com/zyxue/ncbitax2lin) used to add taxonomic information to the local reference database. Requires pip to install.
- gnu-sed (gsed) [Installation Instructions Here](https://formulae.brew.sh/formula/gnu-sed). Nearly all unix-based computing clusters use this as the default (including UI RCDS, FARM and Barbera), but it needs to be installed on a mac and aliased to 'sed'. It can be easily installed with [homebrew](https://brew.sh/)
- R: This pipeline was optimized for R version 4.2.3 
  Packages required:
  - optparse
  - lubridate
  - tidyverse
  - rentrez

**To Install R Packages on the Cluster**
ssh into zaphod, arthur, petunia, whale or another standalone server

`cd ~`

`mkdir Rpackages`

`cd ~/Rpackages`

`module load R/4.2.3`

`R`

When R is running, enter:

`install.packages("optparse",lib="~/Rpackages/",dependencies=T)`

`install.packages("lubridate",lib="~/Rpackages/",dependencies=T)`

`install.packages("tidyverse",lib="~/Rpackages/",dependencies=T)`

`install.packages("rentrez",lib="~/Rpackages/",dependencies=T)`

Or follow [the tutorial here](https://www.hpc.uidaho.edu/compute/Applications/R.html)
    
**Providing an NCBI API Key** 

Currently, you must provide Step 1 with an NCBI API key. Getting an NCBI account is free: [sign up here](https://account.ncbi.nlm.nih.gov/signup/). Once logged in, click your name in the upper right hand corner of the screen, click "Account Settings" and scroll down to the button that says "generate API key". Copy this key into the wrapper script and you'll be good to go.

### Installing the pipeline

**Recommended install on the cluster:**

1. make a new directory called your_project

`mkdir your_project`

2. move all your data files into folders called gene1,gene2...geneN for each gene/primer set you want to analyze. Move these folders into the your_project directory:
   
`mv /gene1/ /path/to/your_project/gene1/`

3. Navigate into your project directory:
   
`cd /path/to/your_project`

4. load the git module
   
`module load git`

5. clone this repository into your folder

`git clone https://github.com/sckieran/Metabarcoding_Pipeline_Waits`

6. move the scripts folder and the steps into the your_project directory
   
`mv ./Metabarcoding_Pipeline_Waits/*.sh .`

`mv  ./Metabarcoding_Pipeline_Waits/scripts/ ./scripts/`


**A note on organizing your data**
To maximize pipeline success, we recommend the following directory structure for your project, although there is flexibility:

<img width="604" alt="Screen Shot 2023-06-30 at 12 19 09 PM" src="https://github.com/sckieran/Metabarcoding_Pipeline_Waits/assets/53580356/00ef8cab-b8d2-45c7-b2f5-19c08fac0cf7">


You should avoid spaces and special characters (^,$,%,@,#,!,*) in the names of your files/folders and check for hidden characters like carriage returns (sometimes displayed as ^M or \r) in your filenames. Taxlist is only required if you're building a ref database. genelist is required. You can name your genes/loci/primer sets anything (again, avoiding spaces and special characters), but must provide a single-line, tab-separated list of the terms in a file. If you're building the local reference database, the program will automatically pull these term from the header of your gene search terms list. Your data files must be separated by gene in folders that correspond to the terms in the "genelist" file. This is even true if you only have one marker.

If you simply provide a project directory (folder) containing the step 1 input files and also containing a directory called gene1 (through ...geneN) that contains your fastas, the pipeline will create the rest of the directories for you within that project directory.

# Usage and Tutorial: 
## SLURM 
You **must** run the wrapper on the SLURM cluster, fortyfive.hpc.uidaho.edu. You must provide the script with your RCDS username (to check when jobs finish) and with a max_jobs value. For RCDS users, the max jobs you can submit is 500. However, running 500 simultaneous jobs, even short ones, can clog up the scheduler. We recommend setting max_jobs to between 50 and 250. The longest step, by far, is the classifier. Read pairing should take ~5 seconds per sample, and collapsing is even shorter. Assigning taxonomy can take ~2-5 seconds per _sequence_, but this filter-free method may produce >500,000 sequences to assign. If you have 500,000 sequences and 250 jobs available, you can expect each jobfile to assign about 2,000 sequences, which should take about an hour. However, you are likely to have some of your jobs held in the queue so you don't consume all of the cluster resources, and can expect the entire process to take ~6 hours, depending on how busy the cluster is. Without this jobfile creation, the process would take about 4 days.

Later functionality will add RRA filters into earlier steps of this process, but are recommended only for experienced users.

## The Wrapper 

The pipeline uses a wrapper script. This makes it easy to run only the steps you're interested in, or to make changes and start somewhere in the middle. However, the pipeline is fairly strict about formatting and filenames, so ensure that you thoroughly understand the naming conventions produced by each step if you plan to run only part of it.

The pipeline is run with `sbatch pipeline_wrapper.sh` and contains parameters you must fill in in order for it to run correctly:

**usage**:
- **dir=**$PWD #this is the containing directory for your project. It should contain your scripts folder, your pipeline_wrapper.sh script, and directories called gene1...geneN with your fastqs. You can use $PWD for this parameter as long as the directory structure is correct.
- **prefix=**your_project #a name for your project, will be used to name your outfiles
- **rlib=**"~/Rpackages" #the location of your R packages from the install process.
- **genelist=**$PWD/genelist #the complete path to your genelist, see example file for format. If your genelist is stored in your project_directory and named genelist, you can use $PWD/genelist for this parameter.
- **taxlist=**$PWD/taxlist #the complete path to your taxlist, see example file for format. If your taxlist is stored in your project_directory and named taxlist, you can use $PWD/taxlist for this parameter.
- **retmax=**20 #How many sequences per taxon to return from NCBI. Default is 20. Recommended values: 5-100. Bigger values are recommended if the taxlist includes genera or higher-order taxa.
- **db_dirr=**reference_database #name (not path) of the directory you want to contain your reference database. Default is 'reference_database'.
- **key=**YOUR_NCBI_KEY #mandatory: your NCBI key for querying NCBI's entrez database. See installation instructions to obtain one. It is free.
- **R1_pattern=**"_R1.fastq" #the pattern that ends your forward read fastq files. Default is "_R1.fastq".
- **R2_pattern=**"_R2.fastq" #the pattern that ends your reverse read fastq files. Default is "_R2.fastq".
- **max_jobs=**200 #the maximum number of SLURM jobs to submit. Default is 200. Max for RCDS users is 500. Play nice with others, even cutting the tasks into only 10 jobs (instead of 1) will reduce processing time to <24 hours even for very large datasets.
- **extra_seqs=**"extra_seqs" #a name for any number of files that include extra sequences not pulled from NCBI that you want included in your database. Full name **must* match the name you give here with this format: **extra_seqs_${gene}_sequences.fasta**. Leave it in your project_directory and it will be automatically included in your reference.
- **filter=**TRUE #do you want to filter your data after BLAST and classification? Default is TRUE.
- **taxa_rra=**0.005 #what should your relative-read abundance cutoff be per taxa? Default is 0.005, or 0.5% of total reads in a sample. See relevant literature in your system/gene for recommended cutoffs.
- **identity_cutoff**=98 #what do you want the percent identity cutoff to be? See relevant literature in your system/gene for recommended cutoffs.
- **minlen=**70 #ASVs shorter than this length are discarded before BLASTing. Set to 1 for no filtering, expect processing times to increase. If you need to set separate minlens for each gene, we recommend pre-trimming in cutadapt (-m $minlen) when you perform your adapter trimming.
- **return_low=**TRUE #should BLAST return all results, or only those that hit above your $percent_identity cutoff? Default is TRUE. If FALSE, sequences with no BLAST hits above threshold are returned as "No Hit".
- **user=**your_username #your RCDS username, for checking if jobs are done.

Once the wrapper is filled out, simply comment out any steps you want to skip (ensuring the mandatory outfiles and directories from those steps exist) and run:

`sbatch pipeline_wrapper.sh`


## Walkthrough of Each Step
### Step One: Fetch FASTAs

**Scripts:**
- step_1_get_seqs_for_database.sh
- query_rentrez.R

**What This StepDoes:**

This step does the following, in order:
1. Parses your arguments and makes relevant output directories
2. Executes the script query_rentrez.R, which:
3. Executess a loop using the `rentrez` R package that searches NCBI for each taxa and each gene term you give it. So the search for the first taxa in the example files would query NCBI for "Aplodontus rufia[ORGN] AND ("12S OR 12s Ribosomal RNA OR 12S RNA OR 12S Mitochondrial"). It does this for each taxa, one a at time. It stores the following:
   * Whether or not there are sequences in genbank for that set of gene terms
   * How many sequences there are
   * The FASTA header and sequence for the first 10 (default) or retmax (set by user) matching sequences in NCBI. These sequences are written to a FASTA file in your out directory with the name "taxa"_"gene"_sequences.fasta
4. Creates a summary file with the above information.

That's it for Step 1. We wanted to make sure that users had the opportunity to add their own, potentially off target samples to their databases For instance, there shouldn't be Homo Sapiens hits to trnL or other plant-specific loci, but they may still show up in your data from contamination in the sequencing lane. Similarly, the European Carp genome (Cyprinus carpio) and lots of COVID sequences are full of Illumina adapters, so including a few can help you identify primer dimers and untrimmed sequences in your data.

### Step Two: Validate FASTAs and Make Ref Database

**Scripts:**

- step_2_make_database.sh
- step_2_p1_rmdups.sh

**What This Script Does:**

This step does the following, in order:
1. Parses arguments and makes out directories
2. Moves your "extra_seqs_gene1_sequences.fasta" files into your "reference_database" directory to be validated and incorporated into your database
3. Loops through each taxon's sequence fastas and queries Entrez to get NCBI-assigned taxids that are associated with each taxon.
4. Adds the NCBI-assigned taxid to the end of the FASTA sequence header (the part that starts with ">") for each sequence. This makes taxonomic assessment easier down the road.
5. Concatenates all taxa FASTAs for each gene into one large FASTA.
6. Checks for (and removes) duplicate sequences with the script step_2_p1_rmdups.sh. This is because sequences with identical headers cause makeblastdb to throw an error
7. Uses the `ncbi-blast` function `makeblastdb` with the following parameters: -db_type nucl (options are "nucl" and "prot") -in your_project_gene1_database_sequences.fasta -out your_project_gene1_reference -parse_seqids -blastdb_version 5
8. Moves all the extra taxa-specific FASTA files into folders to keep things tidy.

End result is a set of reference database files (10 of them, created by makeblastdb) that NCBI blast can use as a reference, along with a fasta file for each gene, all in a folder called reference_database (or a name supplied by you) inside your project folder


### Step Three: Merge forward and reverse reads

**Scripts**

- step_3_pears.sh
- pears.sh

**What This Script Does**

This script does the following, in order:
1. Parses arguments and create out directories for one gene at a time ("gene1" in this tutorial)
3. Grabs all the files in your data directory that match the forward/reverse patterns you provide (default: _R1.fastq and _R2.fastq)
4. Assesses the number of jobs it should create based on your max_jobs and the number of samples in your gene1 directory
5. Makes N lists of samples. Each list contains X samples, where X=(num_samples / max_jobs)
6. Uses the software `pear` to merge your forward and reverse reads.
7. Checks if your jobs are done, and when they are:
8. Outputs the results of `pear`, one for each sample, to sample_paired.assembled.fastq

The end result is a single paired.assembled.fastq file containing merged reads for each sample, along with the reads that were not assembled (_paired.unassembled.forward.fastq and _paired.unassembled.reverse.fastq) and the discarded merged reads (internal pears defaults discard mergers <50 bp, to disable this, or to add additional filtering to pears, edit the file `pears.sh` to include your desired parameters, potential parameters can be found in the pears manual [here](https://cme.h-its.org/exelixis/web/software/pear/doc.html)

### Step Four: Collapse unique sequences into ASVs

**Scripts**

step_4_collapse.sh
run_collapser.sh

**What These Scripts Do**

Step 4 does the following:
1. Creates jobs as in step 3, based on your max_jobs and the number of samples found that match *_paired.assembled.fastq
2. Cleans up your unpaired, unassembled and discarded fastqs
3. Runs fastx-collapser on each sample, which functionally pulls the unique sequences out from each sample as a fasta file, naively and with no error correction, counting the number of reads per ASV.
4. Removes singletons and doubletons. This prevents major ASV overflow issues later.
5. Checks if your jobs are done
6. Cleans up your outfiles, collapsed outfiles have the name sample_collapsed.fasta

### Step Five: Make per-sample sequence files

**Scripts**

step_5_mk_seqfiles.sh

**What This Step Does**

Step 5 does the following, in order:
1. Reformats your collapsed fastas into per-sample ASV tables with the format "sequence  reads", tab-separated
2. Cleans up your outfiles and ensures you're ready for BLAST


### Step Six: BLAST and Classify

**Scripts**

step_6_blast.sh
run_tax.sh

**What This Step Does**

Step 5 does the following, in order:
1. Copies your per-sample sequence files into a directory called "gene1_out"
2. Combines the sequences from all samples into one master file, calls uniques and filters out anything <minlen.
3. Creates a BLAST-formatted fasta from your unique ASV master file. Sequences are given a unique identifier formatted as seq_000001 through seq_999999. File is called your_project_gene1_combined_ASVs.fasta
4. Parses your taxonomy file, or attempts to install and run ncbitax2lin to produce a taxonomy file.
5. BLASTs your BLAST-formatted query fasta against your reference database. If "return_low" is set to TRUE, returns all the hits with the highest bitscore (no matter how many are equally good), plus the hits from the next two highest bitscores. If return_low is set to FALSE, returns only hits with a percent identity above your identity_cutoff. Runs on 4 threads.
6. Parses your raw BLAST results to add convenient taxonomy information (taxid, species) from description of the query sequence.
7. As in steps 3 and 4, assesses how many jobs to create to assign taxonomy based on your # of sequences and creates jobfiles for each job.
For each jobfile:
  8.  Pulls every matching hit from the BLAST output
  9. Assesses it to determine the highest scoring hit(s) (see the previous section for a description of how the program "scores" the hits)
  10. Checks how many species are present in the highest scoring hit(s). If it is only one, the program stores the species and taxonomic information for that sequence as the "best hit". If the top BLAST hit(s) include multiple species, the program then fetches the taxonomy tree for each species among the highest scoring hits and walks up the tree, checking at each level (genus, family, class, order, phylum) whether there is taxonomic agreement among the best hits. If no agreement is found at the phylum level, the sequence is recorded as a no-hit. Otherwise, the lowest commonly-shared taxon is reported as the "best hit" for that sequence.
  11. Returns a "best hit" entry for each sequence that includes the sequence number, the length of the alignment, the percent identity of the top hit, the assigned taxa, the taxonomic information for that hit (phylum, order, class, family, genus), the bitscore of the top hit, the number of species present in the potential best hits, and a comma-separated list of each species that was an equally-good hit.
12. Checks if all jobs are done
13. Concatenates each jobfile's best_hit outfile together, cleans up extra outfiles and gets ready to produce the taxa table.

### Step Seven: Format Unfiltered Taxatable

**Scripts**

step_7_taxatable.sh
run_ttb.sh

**What This Step Does**

1. As in previous steps, creates jobfiles containing X samples per job file where X=$( num_samples / max_jobs)
2. For each sample, assigns taxonomy and "best hit" information to each ASV in the sample_seqs.txt file
3. Checks if jobs are done
4. When jobs are done, concatenates job outfiles to create a master taxatable containing all ASVs for all samples, with taxonomy and best-hit information
5. Cleans up extra outfiles.

### Step eight: Filter Taxatable

**Scripts**

step_8_filter_data.sh
filter_id_taxa.R

**What This Step Does**

1. Feeds your unfiltered taxatable into R
2. Loops over each sample and does the following:
  3. uses the `dplyr` package in R (pre-loaded onto the RCDS R installation, but also included in `tidyverse`) to filter your data by identity_cutoff.
  4. Groups sequences by taxa of best hit. Only sequences with percent identity > identity_cutoff are included in this grouping.
  5. Makes a list of taxa with total reads > cutoff, where cutoff=$( taxa_rra * total_reads_in_sample ), the "passing_taxa" list
  6. Filters the sample's sequence to include sequences from only taxa in the "passing_taxa" list
7. Re-builds the taxatable to include all samples
8. Produces outfile your_project_filtered_taxatable.txt and cleans up outfiles.



