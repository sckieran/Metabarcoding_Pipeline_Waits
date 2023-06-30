# Waits Lab Metabarcoding Training Pipeline
Shannon Blair 2023

A rough pipeline and tutorial to analyze metabarcoding data using dada2. This explains my system. This project is theoretically able to run on the UIdaho RCDS cluster and on a mac running Monterey. It is optimized for cluster use currently. However, it is intended to be an analysis guide for students working on their own projects, rather than a button-press way to receive your output. This is a teaching tool, so it is pretty wordy. Please read through it to understand all the decision points in this analysis.

This pipeline assumes reasonably good resolution of locus data and is designed as a "first pass" at your metabarcoding data, not a finished product. **This pipeline grabs the top BLAST hit for each unique sequence. It does not assess multiple BLAST hits.** There is an optional extension script that will attempt to resolve identically-good BLAST hits. This extension script was built for a specific project where the resolution of the locus was well-understood and there were relatively few possible target taxa. It is extremely naive and basic. If your resolution is low (ie, you have many disparate taxa with perfect or near-perfect matches at your locus), or if you don't understand the resolution of your locus, consider using a Bayesian taxonomic assignment algorithm like the ones found in DADA2 or QIIME instead.

**This pipeline is not designed for microbial metagenomics at 16S**. I say this because there is an immense, incredibly well-maintained array of resources for 16S microbial amplicon sequencing, some commercial and some open source, that will be infinitely better than this pipeline.

## Inputs and Installation

### Inputs
This pipeline requires the following inputs:
- Demultiplexed, appropriately-named paired-end metabarcoding data in fastq (or fastq.gz) format, with adapters trimmed. Each gene or primer set should be analyzed separately, and should be in different folders.
- A parameters file (see template) with your filtering parameters. You don't need this until step 3, so you can fill it out after you've seen the quality report on your reads. It must be in exactly the format as the template.
- The output of a single run of ncbitax2lin. This script will attempt to install and run ncbitax2lin if it can't find a taxonomy file. If you are having issues with permissions and install, you can run it yourself by copying the commands from [this readme](https://github.com/zyxue/ncbitax2lin), the script assumes a name of "ncbi_lineages_[date_of_utcnow].csv.gz" but is agnostic to the date, and you can provide a name if preferred.
- A tab-separated, single-line list of genes/loci you're interested in. If you provide a genelist in step 1, the script will use the header from that list.

Additionally, to run the local database building tool, you need:
- A taxa file containing the scientific names of your target taxa, one per line, with the header "taxname"
- A file with lists of common terms for your target genes. One gene per column, as many permutations on the gene as you'd like, one per line (ie, "Cytochrome Oxidase I", "COI", "COX1"). See the example files for a template. The columns should have a header (which you can repeat in the body) that is ta short, human-readable name of the gene that contains no spaces, slashes, quote marks or other special characters. For example, instead of heading your column "Cytochrome Oxidase I", head it "COI". This file can serve as your "genelist" file for steps 2 and 3 as well.

Data is often demultiplexed by the sequencing service company at no (or minor) cost. However, if your data has not been demultiplexed, we recommend using either [fastq-multx](https://github.com/brwnj/fastq-multx) or the demux-by-name function of [BBMap](https://github.com/BioInfoTools/BBMap). Some demultiplexers (fastq-multx, for instance) do automatic trimming, others do not. The dual-indexed fusion primers often used in metabarcoding may require extra trimming of the overhangs. We recommend using [cutadapt](https://cutadapt.readthedocs.io/en/stable/) to trim. 

### Software and Installation

The pipeline is relatively light on software. R package management will probably be the most intensive thing.

- Command-Line BLAST [Installation Instructions Here](https://www.ncbi.nlm.nih.gov/books/NBK279690/)
- NCBItax2lin: [see here](https://github.com/zyxue/ncbitax2lin) used to add taxonomic information to the local reference database. Requires pip to install.
- gnu-sed (gsed) [Installation Instructions Here](https://formulae.brew.sh/formula/gnu-sed). Most computing clusters use this as the default, but it needs to be installed on a mac. It can be easily installed with [homebrew](https://brew.sh/)
- R: This pipeline was optimized for R version 4.1.2 -- Bird Hippie and has been tested on R version 4.2.3.
- Packages required:
  - optparse
  - tidyverse
  - lubridate
  - rentrez
  - dada2
  - stringr
  - bioconductor

**If you are using the UI RCDS cluster:**
There is a good tutorial for installing r packages here: [link to RCDS](https://www.hpc.uidaho.edu/compute/Applications/R.html)
However, you need to make one change to the tutorial: you must load the newest R module (`module load R/4.2.3`). To install dada2, you must first install bioconductor, [follow the instructions on the dada2 website](https://benjjneb.github.io/dada2/dada-installation.html)

**for this pipeline to work on the RCDS cluster, you must install your packages into a specific place: ~/Rpackages. If you install your R packages somewhere else, you _must_ go into the R scripts and change the `lib="~/Rpackages"` line of each `library` command by hand (or with sed).**

### Before you Begin
The pipeline assumes all internal R scripts are in the same folder as the shell scripts. 

**A note on organizing your data**
To maximize pipeline success, we recommend the following directory structure for your project, although there is flexibility:


                                  |------"sample_1_gene1_R1.fastq"
                                  |
                                  |------"sample_1_gene1_R2.fastq"
                |-------- gene1-|
                |                 |------"sample_2_gene1_R1.fastq"
                |                 |
                |                 |------"sample_2_gene2_R2.fastq"
 Your_Project --                                
                |                  |------"sample_1_gene2_R2.fastq"
                |-------- gene2 -|
                |                  |------"sample_2_gene2_R1.fastq"
                |
                |--------"taxlist"
                |--------"genelist"
                |
                 -------database

You should avoid spaces and special characters (^,$,%,@,#,!,*) in the names of your files/folders and check for hidden characters like carriage returns (sometimes displayed as ^M or \r) in your filenames. Taxlist is only required if you're building a ref database. genelist is required. You can name your genes/loci/primer sets anything (again, avoiding spaces and special characters), but must provide a single-line, tab-separated list of the terms in a file. If you're building the local reference database, the program will automatically pull these term from the header of your gene search terms list. Your data files must be separated by gene in folders that correspond to the terms in the "genelist" file. This is even true if you only have one marker.

## Step 0: Build a Local Reference Database

We recommend a local reference database to improve the accuracy of detections. Future versions of this pipeline may also include a method to format this library so it can be used with DADA2's taxonomic assignment algorithm. This script uses [rentrez](https://cran.r-project.org/web/packages/rentrez/vignettes/rentrez_tutorial.html), an excellent R package for querying NCBI. However, you can skip this step and use remote BLAST instead. Building a reference library may take several hours, but the actual searches are much faster than remote BLAST.

## Things to Consider
- This script assumes a relatively small reference library, <500 taxa at <5 genes. This script is rate limited to avoid overloading the NCBI API, but can still run into odd errors when internet connections are unstable, and doesn't always play nicely with VPNs. If you're getting weird errors like "unexpected EOF" or "http 400", try disconnecting any VPNs or checking firewalls. I am generally unable to help you troubleshoot these issues.
- We recommend that you include common contaminants in your database for better error detection. At minimum, we recommend *Homo sapiens* for any genes that might amplify vertebrates and *Cyprinus carpio* for all projects. As of mid-2023, the *Cyprinus carpio* genome is full of Illumina adapters, so if you have a lot of primer dimer in your reads, or if you fail to trim your adapters, you will have many hits to this species. I assume it will get a new genome at some point, so this may not be true forever. 
- This script can search at any taxonomic level. For species, include the binomial name. For subspecies, check NCBI for the correct name of the subspecies you're interested in, or use taxids. For genera and above, use the single-word taxon name.
- There is an extension to the R script that can check for sequences within a genus for any species that has no available relevant sequences in genbank, to add congeners and therefore accuracy to the database. By default, it is commented out. To use it, simply remove the comment character (#) inside the query_rentrez.R script, starting at the section labeled "search no-hits by genus" and going to the end of the file. **note: currently, this will not work.**
- This script doesn't check for the wrong organism in the search results. This allows the script to search for organisms at any taxonomic rank. However, there is a version of this script called step_0_make_local_database_by_taxid.sh that you can use instead if you are worried about sister taxa slipping into your search results. The input for that script is a taxfile containing a list of taxids, which you are responsible for supplying. In theory, you can use [this NCBI tool](https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi) to look up taxids. **The taxid version of this script only accepts species/subspecies as input, not genera or any higher taxonomic rank, and will check to ensure that only the target taxid is included in the output files.** **note: the taxid script isn't ready yet**
- 
### Inputs
- A taxa file containing the scientific names of your target taxa, one per line, with the header "taxname"
- A tab-delimited file with lists of genes with "gene terms". Each column corresponds to a gene, each row is different ways to describe that gene. This just provides the best search radius. So 12S would include "12S" "12S Ribosomal RNA" "12S Mitochondrial Sequence". The columns can be of different lengths. We recommend 2-3 terms per gene for the most common metabarcoding loci (12S, COI, ITS, 16S)


### Usage
**Cluster Usage**
To use this script on the RCDS 44 cluster, edit the step_1_wrapper.sh file to include your relevant project names and paths. To specify your current working directory  as your project directory, you can use $PWD. For example, in the wrapper, `dirr=$PWD` would tell the script to use your current directory as the project directory. Similarly, `genelist=$PWD/genelist` will tell the script to look for a file called "genelist" in your current directory.

Then run
`sbatch step_3_wrapper.sh`

Make sure the step_1_wrapper, and the step_1.sh script are both in your directory.

**to run on a head node**

**Arguments/Options**
The script requires four command line arguments:
* -n a name for the output files, ex. "your_project"
* -t the path to the file containing a list of taxa
* -g the path to the file containing the gene terms
* -d the path to the output directory

**Optional Arguments:**
* -r (retmax) maximum number of search records to return. Default 10. Setting this value very high (>50) may cause problems with the rentrez search and will probably add off-target sequences (wrong gene more likely than wrong organism) to your database. I recommend getting an NCBI API key for large databases, see [the rentrez tutorial](https://cran.r-project.org/web/packages/rentrez/vignettes/rentrez_tutorial.html#rate-limiting-and-api-keys) for more info.
* -c combine: choice of "comb" or "sep" (default). If "comb" is invoked, the fastas from all genes will be combined into a single database. This is not recommended. Sequences from different genes/loci must be processed separately (because they will almost certainly be different lengths)
  
`bash step_0_make_local_database.sh -n your_name -t /path/to/taxfile -g /path/to/genefile -d /path/to/outputdir -r 10 -c sep -r ~/Rpackages`

### Output

In your output directory:
- A tab-delimited summary file with information about each taxon including the taxid, number of available sequences on NCBI (at any locus), whether sequences are available for any of your target loci, and a yes/no for each locus, called "your_project_database_taxa_summary.txt"
- A list of taxa that had no hits for your relevant genes in NCBI, called "your_project_database_taxa_no_hits.txt"
- A fasta file for each gene/locus called your_project_gene_database_seqs.fasta
- Directories for each gene containing the fasta files for each taxa separately, which can be deleted if you're done with them. 
- A searchable local NCBI reference database for each gene. These consists of 10 files, each with the name "your_project_gene_reference.[suffix]".


**If you already have a reference library:**
Your reference library needs to be an NCBI blast database (made with makeblastdb). You can use your own fasta file to create it and name it whatever you'd like. However, please note that only reference libraries made using NCBI accession numbers and the -parse_seqids option of makeblastdb will be able to take advantage of the part of this pipeline that assessess equally-good BLAST hits to find a consensus taxon.

If you only have a few extra sequences with no accession numbers (for instance, that you sequenced yourself), there is a somewhat-tedious workaround for this. The fasta header line for that sequence should be edited to have the following format: >unique_identifier species name taxid=unique_value. It should be space-delimited. Then, to the bottom of your NCBItax2lin output file, you can add a row with your unique_value and the relevant taxonomic information. Your unique value and unique_identifier can be anything that doesn't contain spaces or special characters. The format for the taxonomy info in the NCBItax2lin output is visible in the header, and is comma-separated.

## Step 2: Data Quality Analysis

### Inputs
- Your demultiplexed, appropriately-named, trimmed, paired-end data files. Each gene should have data files in a separate folder named "gene1" "gene2" etc for each of your genes/loci/primer sets.
**What does this mean??**
  -  **demultiplexed** means that your reads have been separated into separate files for each sample.
  -  **appropriately-named** means that your files are named with the following scheme ${name}${pattern}. For example, "Sample-1_R1.fastq". Here the name would be "Sample-1" and the pattern would be "_R1.fastq". You would then tell the program that your _pattern1_ is "_R1.fastq". Note: When ${pattern} is stripped from your filenames, **all filenames in your data folder must be unique**. DADA2 cannot process duplicated sample names and will assume that sample names are whatever is left when ${pattern} is stripped off the filename, so make sure you don't have any naming anomalies in your data. Pay special attention to your positives and negatives, which is sometimes where these errors occur. The defaults for ${pattern1} and ${pattern2} are "_R1.fastq" and "_R2.fastq".
  -  **paired-end** As of right now, this pipeline can only handle paired-end data. This means you should have a forward and a reverse read file for each sample.
-  a "genelist" file, a single-row, tab separated file with gene/locus names that match the containing folders for your data files. If you did local database building, this can be the same genelist you used there.

### Usage

**Running on the cluster**
Fill out the step_2_wrapper.sh file just like you did for step 1, then run:
`sbatch step_2_wrapper.sh`

Once it's done running, adjust your filtering_params file accordingly. You're ready for step 3.
