# Waits Lab Metabarcoding Training Pipeline
Shannon Blair 2023

A rough pipeline and tutorial to analyze metabarcoding data using dada2 with a limited local reference database. It is intended to be an analysis guide for students working on their own projects.

This pipeline assumes reasonably good resolution of locus data and is designed as a "first pass" at your metabarcoding data, not a finished product. The way this pipeline assesses BLAST hits is both **naive and conservative**. If more than one equally-likely BLAST hit is available (judged by percent identity), it will attempt to assign taxonomy to each hit using the output of ncbitax2lin, then walk up the taxonomy tree until it finds a consensus rank. If no consensus is available at the phylum level, it calls a "no-hit". Therefore **this pipeline only works with reference libraries made using this pipeline.** However, you can use your own FASTAs as input to step 2 of this pipeline if you want to skip querying rentrez, as long as each sequence in that FASTA begins with a header line in the format >[Valid NCBI Accession Number]. 

**What if I want my database to include every organism available for X gene?**
This pipeline is for limited-taxa reference databases or remote querying of NCBI. It is not designed to manage the curation of a large database that includes, say, all inverts with COI sequences in NCBI. If you want a larger database, consider using another reference database building tool, for instance [RESCRIPt](https://github.com/bokulich-lab/RESCRIPt) or [bcDatabaser](https://bcdatabaser.molecular.eco/). This pipeline is being actively developed to better improve performance with other databases. Any FASTA file that contains unique sequences each with a valid NCBI accession will work in step 2 of the pipeline, but runtimes for step 2 will be very long for huge databases (>1000 taxa).

**What if I have a few of my own in-house sequences I want to add to my database?**
That's fine, but it will add some tedioius work upfront. First, check if the species you are including have taxIDS assigned in NCBI. All taxa with at least one sequence in the NCBI nucleotide or sra database have a taxid, and you can look it up here: [here](https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi). Put your homebrewed sequences in a file called extras_gene1_sequences.fasta with the header format ">unique_identifier species name taxid=unique_taxid" and put it in your reference_database directory. If your species has a taxID assigned in NCBI, use that as the unique_taxid. Otherwise, you can assign one yourself, I recommend a series of letters and numbers separated by underscores. Do not use a strictly numeric value, you will almost definitely accidentally assign a real taxID to your species. You can assign any unique_identifier to the sequence (to replace the accession number), but it must not include spaces or special characters and each sequence must be uniquely named. Next, also in your reference database directory, install and run ncbitax2lin. Then, unzip and modify the ncbitax2lin file. Add the unique_taxids you assigned to the end of the file and add taxonomic information for each species following the comma-separated format of the ncbitax2lin file (described in the header of the file). This is tedious to do by hand, so we really recommend submitting your homebrew sequences to NCBI if they're long enough!

**What if I want to use remote BLAST to query all of NCBI?**
Sure, that's fine, start the pipeline at step 3 and run the step_4_remote_blast_option.sh or step_4_remote_wrapper.sh (for head nodes and SLURM clusters, respectively). Expect the BLAST search to take several hours depending on the number of sequences, and make sure you have a good internet connection.

**This pipeline is not designed for microbial metagenomics at 16S**. I say this because there is an immense, incredibly well-maintained array of resources for 16S microbial amplicon sequencing, some commercial and some open source, that will be infinitely better than this pipeline.

## Inputs and Preparation

### Inputs
This pipeline requires the following inputs:
- Demultiplexed, appropriately-named paired-end metabarcoding data in fastq (or fastq.gz) format, with adapters trimmed. Each gene or primer set should be analyzed separately, and should be in a different folder called gene1, gene2...geneN.
- A parameters file (see template) with your filtering parameters. You don't need this until step 3, so you can fill it out after you've seen the quality report on your reads. It must be in exactly the format as the template.
- The output of a single run of ncbitax2lin. This script will attempt to install and run ncbitax2lin if it can't find a taxonomy file. If you are having issues with permissions and install, you can run it yourself by copying the commands from [this readme](https://github.com/zyxue/ncbitax2lin), the script assumes a name of "ncbi_lineages_[date_of_utcnow].csv.gz" but is agnostic to the date, and you can provide a name if preferred.
- A tab-separated, single-line list of genes/loci you're interested in. If you provide a genelist in step 1, the script will use the header from that list.

Additionally, to run the local database building tool, you need:
- A taxa file containing the scientific names of your target taxa, one per line, with the header "taxname"
- A file with lists of common terms for your target genes. One gene per column, as many permutations on the gene as you'd like, one per line (ie, "Cytochrome Oxidase I", "COI", "COX1"). See the example files for a template. The columns should have a header (which you can repeat in the body) that is ta short, human-readable name of the gene that contains no spaces, slashes, quote marks or other special characters. For example, instead of heading your column "Cytochrome Oxidase I", head it "COI". This file can serve as your "genelist" file for steps 2 and 3 as well.

Data is often demultiplexed by the sequencing service company at no (or minor) cost. However, if your data has not been demultiplexed, we recommend using either [fastq-multx](https://github.com/brwnj/fastq-multx) or the demux-by-name function of [BBMap](https://github.com/BioInfoTools/BBMap). Some demultiplexers (fastq-multx, for instance) do automatic trimming, others do not. The dual-indexed fusion primers often used in metabarcoding may require extra trimming of the overhangs. We recommend using [cutadapt](https://cutadapt.readthedocs.io/en/stable/) to trim. 

### Software and Installation

You will need
- Command-line github (pre-installed on UI RCDS Cluster, UCD FARM and Barbera)
- Command-Line BLAST [Installation Instructions Here](https://www.ncbi.nlm.nih.gov/books/NBK279690/) (pre-installed on UI RCDS Cluster, UCD FARM and Barbera)
- NCBItax2lin: [see here](https://github.com/zyxue/ncbitax2lin) used to add taxonomic information to the local reference database. Requires pip to install.
- gnu-sed (gsed) [Installation Instructions Here](https://formulae.brew.sh/formula/gnu-sed). Nearly all unix-based computing clusters use this as the default (including UI RCDS, FARM and Barbera), but it needs to be installed on a mac and aliased to 'sed'. It can be easily installed with [homebrew](https://brew.sh/)
- R: This pipeline was optimized for R version 4.2.3 
- Packages required:
  - optparse
  - tidyverse
  - lubridate
  - rentrez
  - stringr
  - bioconductor
  - dada2 (see dadad2 installation instructions [here](https://benjjneb.github.io/dada2/dada-installation.html) for bioconductor installation of dada2, may need to remove the "version=" flag to install)

**If you are using the UI RCDS cluster:**
There is a good tutorial for installing r packages here: [link to RCDS](https://www.hpc.uidaho.edu/compute/Applications/R.html)
However, you need to make one change to the tutorial: you must load the newest R module (`module load R/4.2.3`). 

**for this pipeline to work on the RCDS cluster, you must install your packages into a specific place: ~/Rpackages. If you install your R packages somewhere else, you _must_ go into the R scripts and change the `lib="~/Rpackages"` line of each `library` command by hand (or with sed).**

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
`mv ./Metabarcoding_Pipeline_Waits/*.sh .
`mv  ./Metabarcoding_Pipeline_Waits/scripts/ ./scripts/`


**A note on organizing your data**
To maximize pipeline success, we recommend the following directory structure for your project, although there is flexibility:

<img width="604" alt="Screen Shot 2023-06-30 at 12 19 09 PM" src="https://github.com/sckieran/Metabarcoding_Pipeline_Waits/assets/53580356/00ef8cab-b8d2-45c7-b2f5-19c08fac0cf7">


You should avoid spaces and special characters (^,$,%,@,#,!,*) in the names of your files/folders and check for hidden characters like carriage returns (sometimes displayed as ^M or \r) in your filenames. Taxlist is only required if you're building a ref database. genelist is required. You can name your genes/loci/primer sets anything (again, avoiding spaces and special characters), but must provide a single-line, tab-separated list of the terms in a file. If you're building the local reference database, the program will automatically pull these term from the header of your gene search terms list. Your data files must be separated by gene in folders that correspond to the terms in the "genelist" file. This is even true if you only have one marker.

## Optional: Step 1: Build a Local Reference FASTA File

We recommend a local reference database to improve the accuracy of detections. Future versions of this pipeline may also include a method to format this library so it can be used with DADA2's taxonomic assignment algorithm. This script uses [rentrez](https://cran.r-project.org/web/packages/rentrez/vignettes/rentrez_tutorial.html), an excellent R package for querying NCBI. 

### Things to Consider
- This script assumes a relatively small reference library, <500 taxa at <5 genes. This script is rate limited to avoid overloading the NCBI API, but can still run into odd errors when internet connections are unstable, and doesn't always play nicely with VPNs. If you're getting weird errors like "unexpected EOF" or "http 400", try disconnecting any VPNs or checking firewalls. I am generally unable to help you troubleshoot these issues.
- We recommend that you include common contaminants in your database for better error detection. At minimum, we recommend *Homo sapiens* for any genes that might amplify vertebrates and *Cyprinus carpio* for all projects. As of mid-2023, the *Cyprinus carpio* genome is full of Illumina adapters, so if you have a lot of primer dimer in your reads, or if you fail to trim your adapters, you will have many hits to this species. I assume it will get a new genome at some point, so this may not be true forever. 
- This script can search at any taxonomic level. For species, include the binomial name. For subspecies, check NCBI for the correct name of the subspecies you're interested in, or use taxids. For genera and above, use the single-word taxon name.

- 
### Inputs
- A taxa file containing the scientific names of your target taxa, one per line, with the header "taxname"
- A tab-delimited file with lists of genes with "gene terms". Each column corresponds to a gene, each row is different ways to describe that gene. This just provides the best search radius. So 12S would include "12S" "12S Ribosomal RNA" "12S Mitochondrial Sequence". The columns can be of different lengths. We recommend 2-3 terms per gene for the most common metabarcoding loci (12S, COI, ITS, 16S)


### Usage
**Cluster Usage**
To use this script on the RCDS 44 cluster, edit the step_1_wrapper.sh file to include your relevant project names and paths. To specify your current working directory  as your project directory, you can use $PWD. For example, in the wrapper, `dirr=$PWD` would tell the script to use your current directory as the project directory. Similarly, `genelist=$PWD/genelist` will tell the script to look for a file called "genelist" in your current directory.

Then run
`sbatch step_1_wrapper.sh`

Make sure the step_1_wrapper, and the step_1.sh script are both in your directory.

**to run on a head node**

**Arguments/Options**
The script requires four command line arguments:
* -n a name for the output files, ex. "your_project"
* -t the path to the file containing a list of taxa
* -g the path to the file containing the gene terms
* -d the path to your project directory
* -h the **name** (not path) you want to give to your output directory. Path with be /project_directory/output_directory. Default is "reference_database".

**Optional Arguments:**
* -r (retmax) maximum number of search records to return. Default 10. Setting this value very high (>50) may cause problems with the rentrez search and will probably add off-target sequences (wrong gene more likely than wrong organism) to your database. I recommend getting an NCBI API key for large databases, see [the rentrez tutorial](https://cran.r-project.org/web/packages/rentrez/vignettes/rentrez_tutorial.html#rate-limiting-and-api-keys) for more info.
  
`bash step_0_make_local_database.sh -n your_name -t /path/to/taxfile -g /path/to/genefile -d /path/to/outputdir -r 10 -c sep -r ~/Rpackages`

### Output

In your output directory:
- A tab-delimited summary file with information about each taxon including the taxid, number of available sequences on NCBI (at any locus), whether sequences are available for any of your target loci, and a yes/no for each locus, called "your_project_database_taxa_summary.txt"
- A list of taxa that had no hits for your relevant genes in NCBI, called "your_project_database_taxa_no_hits.txt"
- A fasta file for each gene/locus called your_project_gene_database_seqs.fasta


**If you already have a reference library:**
Your reference library needs to be an NCBI blast database (made with makeblastdb) and we **strongly recommend** re-building your reference library using this pipeline in order to avoid errors downstream. If you wish to build your reference library using a FASTA you made yourself, or if you wish to add extra sequences to your library (such as off-target taxa with no representation at your gene), simply name the fasta any_name_gene1_sequences.fasta. You can have as many fastas as you want, just put them in the database directory you want (or the one you specified in step 1). The second step of this pipeline will add taxIDs to your sequences as long as each sequence begins >[Valid NCBI Accession]. You must have a copy of each fasta that you want included for each gene1...geneN that you're analyzing. 

## Step 2: Build the local reference database
This script will grab everything in your database directory called *_gene1_sequences.fasta for each gene1...geneN, add taxIDs based on the NCBI accession number, check for duplicate entries, and then make a blast database. As noted above, you can supply your own FASTA and start the pipeline here, but your FASTA MUST use real NCBI accession #s in order for the taxid lookup to work. Theoretically, the pipeline should still work when a non-standard FASTA is supplied, just without the taxonomic decisionmaking, but behavior may be unpredictable.

### Inputs
- A tab-delimited genelist file (see step 1)
- The FASTA files output from step 1, plus any other fastas you want to add, with the naming format *_gene1_sequences.fasta for gene1...geneN

### Usage

**Cluster usage**
Fill in the blanks in the step_2_wrapper.sh script, including the path to your project directory and the **name** (not path) of the reference_database directory you chose. This directory must be nested inside your project directory, and it must contain at least 1 fasta file with the name scheme *_gene1_sequences.fasta for gene1...geneN.

dirr=/path/to/project/directory (can be $PWD if you're running the script from inside the project directory)
db_dirr=reference_database [name, not path, of database directory]
genelist=/path/to/genelist 
taxlist=/path/to/taxlist
prefix="your_project"
retmax=10 [optional: leave blank for default. Default=10].

`sbatch step_1_wrapper.sh`

**To Run on the Head Nodes**
**Arguments/Options**
The script requires four command line arguments:
* -n a name for the output files, ex. "your_project", should be the same name as step 1.
* -g the path to the file containing the gene terms
* -d the path to the project directory
* -h the **name** (not path) to the output (reference database) directory you chose in step 1. Default is "reference_database".
* -r 10 the number of matches you want to return. More = a longer runtime but potentially higher diversity coverage in your FASTA files. We recommend <100 to help manage NCBI API limits.
  
### Outputs
- Directories for each gene containing the fasta files for each taxa separately, which can be deleted if you're done with them. 
- A searchable local NCBI reference database for each gene. These consist of 10 files, each with the name "your_project_gene_reference.[suffix]".

## Step 3: Data Quality Assessment

### Inputs
- Your demultiplexed, appropriately-named, trimmed, paired-end data files. Each gene should have data files in a separate folder named "gene1" "gene2" etc for each of your genes/loci/primer sets.
**What does this mean??**
  -  **demultiplexed** means that your reads have been separated into separate files for each sample.
  -  **appropriately-named** means that your files are named with the following scheme ${name}${pattern}. For example, "Sample-1_R1.fastq". Here the name would be "Sample-1" and the pattern would be "_R1.fastq". You would then tell the program that your _pattern1_ is "_R1.fastq". Note: When ${pattern} is stripped from your filenames, **all filenames in your data folder must be unique**. DADA2 cannot process duplicated sample names and will assume that sample names are whatever is left when ${pattern} is stripped off the filename, so make sure you don't have any naming anomalies in your data. Pay special attention to your positives and negatives, which is sometimes where these errors occur. The defaults for ${pattern1} and ${pattern2} are "_R1.fastq" and "_R2.fastq". If your data are gzipped, you'll need to change this.
  -  **paired-end** As of right now, this pipeline can only handle paired-end data. This means you should have a forward and a reverse read file for each sample.
-  a "genelist" file, a single-row, tab separated file with gene/locus names that match the containing folders for your data files. If you did local database building, this can be the same genelist you used there.

### Usage

**Running on the cluster**
Fill out the step_3_wrapper.sh file:

`dirr=$PWD #path to your project directory. If you're running this code from inside your directory (recommended), you can set this to $PWD.
genelist=$PWD/genelist #path to your genelist, a tab-separated single-line list of genes/primer sets that correspond to your data folders.
prefix="your_project" #for ease of use, this should match the name you gave in Step 1.
pattern1="_R1.fastq"
pattern2="_R2.fastq"
num_graphs=24 #change this based on how many quality profiles you want to look at. Adding more adds computational time to this script.`


then run:

`sbatch step_3_wrapper.sh`

**Running in the head node**
The  program takes the following arguments:
* -n a name for your outfiles. Should match step 1 for easiest working.
* -g the path to the list containing your gene terms.
* -d the directory that contains your project, probbaly the directory you're running the code from ($PWD)
* -p the pattern for the end of your R1 files. Default is "_R1.fastq".
* -q the pattern for th ened of your R2 files. Default is "_R2.fastq".
* -k the number of samples you want to view quality plots for. Samples will be randomly selected . Default is 24. Increasing this value increaeses computational time.
* -l path to your R packages library

`bash step_3_quality_check_reads.sh	-n name -g genelist -d directory -p pattern1 -q pattern2 -k num_graphs -l rlib`

### Outputs

Step 3 will output two PDFs into your_directory/reports/. They will contain k quality plots from dada's plot_quality_profile function. 

Once it's done running, evaluate the read quality. Look at where the quality drops off the read and consider your amplicon and overlap sizes. Then, adjust your params_file accordingly. You're ready for step 4.

## Step 4: Filtering data, performing ASV selection, blasting ASVs and building taxtables.

Step 4 is a big step. First it will filter your data using the parameters you specify in the params_file. You should have one params_file for each gene (as you should filter your different loci each according to their length and sequencing quality). Your files should be named "params_file_gene1" and "params_file_gene2" for each gene. You can specify the name of "params_file" in the wrapper/command line, as long as the name is followed by "_geneN" for each gene.

Most of the filtering parameters are optional and can be left blank to use dada2 defaults. However, the project directory, relative read abundance cutoff (set to 0 for no RRA filtering), filter directory and multithread options must have some input. Multithread must be TRUE or FALSE (set to TRUE for clusters and macs). **current behavior ignores user-input filtering directory and puts all outs into a directory called "filtered".**

After filtering, the script will output _seqs.txt files that contain the unique ASVs retained for each sample, and their read counts. It will also output a raw ASV table (ASVs as rows and samples as columns, read counts as data) in a directory called `/results_tables/`. Then it will create a fasta file with each unique ASV present in any of your samples. Then, it blasts those sequences against the reference database you made in step 1 and does some calculations about your hit rate. A very low hit rate suggests you should consider filtering your data more stringently, or expanding your reference database.

After that, it will look for a taxonomy file. If you have one already downloaded from NCBI, it will skip trying to install and download. Otherwise, it will attempt to install and download. If you run into permissions issues, simply follow the steps to run ncbi2tax yourself in your database folder and then run step 3 again. It will find it and skip the install step. 

After adding taxonomy to all hits, it will evaluate equally-identical blast hits and pick the taxonomic rank at which there is agreement. This is a naive (but conservative) way of assessing hits. The raw_blast_out and raw_blast_out_with_tax files will be available in the out directory for you to examine by hand.

It then formats and outputs a taxa table, which includes all ASVs for all samples and reports reads, identity, taxa of best hit, and higher taxonomic information. This is a great table to take into excel to look at pivottables, visualize graphs, etc.

### Inputs
- Your sequence data files in folders labeled 'gene1'...'geneN'
- A params_file for each gene, called params_file_gene1...params_file_geneN. You can specify the 'params_file' part of the name.

### Usage

**Command line Usage**
Fill in the wrapper as with step 3. The name and directory must be the same between the steps.

`dir=$PWD #path to the your_project directory. If you are running this code from inside that directory (recommended), you can put $PWD here.
prefix=your_project
db_dir=reference_database #name, not path, of database directory (path will be ${dir}/${db_dir}). Default is 'reference_database'
pattern1="_R1.fastq" #default is _R1.fastq, leave blank for default
pattern2="_R2.fastq" #default is _R2.fastq, leave blank for default
genelist=$PWD/genelist #path to genelist
params_file="params_file" #name of params_file prefix. All parameter filenames should start with this prefix and end with _gene1,_gene2...geneN for each gene/primer set in the genelist. Default is params_file
localdat= #only put something here if your local database is not named "yourproject_gene_reference", which is the default for steps 1 and 2.`

`sbatch step_4_wrapper.sh`

**Head Node Usage**

The program takes the following arguments:
* -n  'your_project' mandatory:  must be the same as your step 2 name
* -g  'genelist' mandatory: same genelist as steps 1 and 2
* -d  '/path/to/your/project/directory' mandatory: should be the same as steps 1 and 2
* -m 'params_file' mandatory: just the first part of the file name, program assumes that it ends '_gene1'...'geneN' for each gene in the genelist. Default is "params_file"
* -p 'pattern1' default: same as in step 2, default is "_R1.fastq"
* -q 'pattern2' default: same as in step 2, default is "_R2.fastq"
* -r name, not path to, folder that contains your reference database (and your NCBI2tax run, if you did that separately). Default is "reference_database"
* -l path to your R packages folder

`bash step_4_filter_and_blast_local.sh -n name -g genelist -d project_directory -m params_file_prefix -p pattern1 -q pattern2 -r database_directory -b ref_database_name -l rlib`

### Outputs

- A folder in your project_directory called "gene1_dada_out" for gene1...geneN. Within that:
    - A folder called "subsampled" containing your subsampled reads for improving error rate estimation time.
    - A folder called "filtered" containing your filtered reads. Altered based on your input to the params file
    - A folder called "sample_seqfiles" containing your _seqs.txt file of unique ASVs and read counts for each sample
    - your_project_gene1_combined_ASVs.fasta, a fasta containing every unique ASV that passed the dada2 filters for any sample, this is the input to blast. Sequences are named seq_1...seq_n.
    - your_project_gene_raw_blast_out, the raw output of the blastn. Contains every hit, the seq ID, the subject accession number, identity and length. Does not report no hits, which are inferred at a later step.
    - your_project_gene_raw_blast_out_with_tax, the raw blast output that adds taxonomy to each blast hit
    - your_project_best_blast_hits.txt a file that has one row for each input sequence (seq_1...seq_n), reporting the best blast hit based on the taxonomic assesment done by this script.
    - your ncbi taxonomy files
- A folder in your project_directory called "results_tables", which includes:
  - a taxa table with every unique ASV reported for every sample, with read count, identity, and taxonomic information
  - a raw ASV table that has each unique ASV (unfiltered for RRA) and read count per sample. Pre-blast, so no taxa information.
- A file in the "reports" folder that tracks the different filtering steps taken by DADA2 and the reads per sample at each step. Only includes samples that had >0 reads pass the initial quality filter.
- 
Yay! You're done! Time to analyze your data in R or excel, affiliate sample names with metadata, and evaluating your parameters. Hooray!

## Alternate Step 4: Filtering data, performing ASV selection, and performing REMOTE BLAST
Very similar to regular step 4, but ignores any local database and performs remote blast. Remote BLAST is very slow for large #s of sequences, so be prepared to wait. This script is currently being tested, and may be revised to make it easier to BLAST large numbers of sequences.

### Inputs
- Your sequence data files in folders labeled 'gene1'...'geneN'
- A params_file for each gene, called params_file_gene1...params_file_geneN. You can specify the 'params_file' part of the name.

### Usage

**Command line Usage**
Fill in the wrapper as with step 4 local. The name and directory must be the same between the steps. If you've done an NCBItax2lin run, put the output in your project directory.

`sbatch step_4_remote_wrapper.sh`

**Head Node Usage**

The program takes the following arguments:
* -n  'your_project' mandatory:  must be the same as your step 2 name
* -g  'genelist' mandatory: same genelist as steps 1 and 2
* -d  '/path/to/your/project/directory' mandatory: should be the same as steps 1 and 2
* -m 'params_file' mandatory: just the first part of the file name, program assumes that it ends '_gene1'...'geneN' for each gene in the genelist. Default is "params_file"
* -p 'pattern1' default: same as in step 2, default is "_R1.fastq"
* -q 'pattern2' default: same as in step 2, default is "_R2.fastq"
* -l path to your R package library

`bash step_3_filter_and_blast.sh -n name -g genelist -d project_directory -m params_file_prefix -p pattern1 -q pattern2 -l rlib

### Outputs

- A folder in your project_directory called "gene1_dada_out" for gene1...geneN. Within that:
    - A folder called "subsampled" containing your subsampled reads for improving error rate estimation time.
    - A folder called "filtered" containing your filtered reads. Altered based on your input to the params file
    - A folder called "sample_seqfiles" containing your _seqs.txt file of unique ASVs and read counts for each sample
    - your_project_gene1_combined_ASVs.fasta, a fasta containing every unique ASV that passed the dada2 filters for any sample, this is the input to blast. Sequences are named seq_1...seq_n.
    - your_project_gene_remote_blast_out, the raw output of the blastn. Contains every hit, the seq ID, the subject accession number, identity and length. Does not report no hits, which are inferred at a later step.
    - your_project_gene_remote_blast_out_with_tax, the raw blast output that adds taxonomy to each blast hit
    - your_project_remote_best_blast_hits.txt a file that has one row for each input sequence (seq_1...seq_n), reporting the best blast hit based on the taxonomic assesment done by this script.
    - your ncbi taxonomy files
- A folder in your project_directory called "results_tables", which includes:
  - a taxa table with every unique ASV reported for every sample, with read count, identity, and taxonomic information
  - a raw ASV table that has each unique ASV (unfiltered for RRA) and read count per sample. Pre-blast, so no taxa information.
- A file in the "reports" folder that tracks the different filtering steps taken by DADA2 and the reads per sample at each step. Only includes samples that had >0 reads pass the initial quality filter.


## Optional: Step 5: Re-run low % hits with NCBI remote BLAST

If your local database is too small or doesn't contain enough contaminant breadth, it might be useful to remotely blast your worst hits (those with low percent identity). That will help you understand if there's more data you can extract from this dataset. This optional step takes an identity percent threshold cutoff you provide, extracts sequences whose **best** BLAST hit in your local database was lower than your cutoff (not inclusive), and re-runs those samples with local blast.

The outputs are a new taxatable with only the remote hits, a new remote_best_hits.txt file with your best hits, and a file called local_vs_remote_best_hits.txt which compares the identity and taxa of the best hits for each sequence with a best local hit below your cutoff.



