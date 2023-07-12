# Waits Lab Metabarcoding Training Pipeline
Shannon Blair 2023

## Overview of Pipeline

A rough pipeline and tutorial to analyze metabarcoding data using dada2 with a limited local reference database. It is intended to be an analysis guide for students working on their own projects.


This pipeline assumes reasonably good resolution of locus data and is designed as a "first pass" at your metabarcoding data, not a finished product. The way this pipeline assesses BLAST hits is both **naive and conservative**. If more than one equally-likely BLAST hit is available (judged by percent identity), it will attempt to assign taxonomy to each hit using the output of ncbitax2lin, then walk up the taxonomy tree until it finds a consensus rank. If no consensus is available at the phylum level, it calls a "no-hit". Therefore **this pipeline only works with reference libraries made using this pipeline.** However, you can use your own FASTAs as input to step 2 of this pipeline if you want to skip querying rentrez, as long as each sequence in that FASTA begins with a header line in the format >[Valid NCBI Accession Number]. 

### Step One: Fetch FASTAs

**Scripts:**
step_1_wrapper.sh
- step_1_get_seqs_for_database.sh
- query_rentrez.R

**What This Script Does:**
This script does the following, in order:
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
step_2_wrapper.sh
- step_2_make_database.sh
- step_2_p1_rmdups.sh

**What This Script Does:**
This script does the following, in order:
1. Parses arguments and makes out directories
2. Loops through each taxa's sequence fastas and queries Entrez to get NCBI-assigned taxids that are associated with each taxon.
3. Adds the NCBI-assigned taxid to the end of the FASTA sequence header (the part that starts with ">") for each sequence. This makes taxonomic assessment easier down the road.
4. Concatenates all taxa FASTAs for each gene into one large FASTA.
5. Checks for (and removes) duplicate sequences with the script step_2_p1_rmdups.sh. This is because sequences with identical headers cause makeblastdb to throw an error
6. Uses the `ncbi-blast` function `makeblastdb` with the following parameters: -db_type nucl (options are "nucl" and "prot") -in your_project_gene1_database_sequences.fasta -out your_project_gene1_reference -parse_seqids -blastdb_version 5
7. Moves all the extra taxa-specific FASTA files into folders to keep things tidy.

End result is a set of reference database files (10 of them, created by makeblastdb) that NCBI blast can use as a reference, along with a fasta file for each gene, all in a folder called reference_database (or a name supplied by you) inside your project folder


### Step Three: Check the quality of your reads

**Scripts**
step_3_wrapper.sh
- step_3_quality_check_reads.sh
- quality_check_reads_in_DADA2.R

**What This Script Does**
This script does the following, in order:
1. Parses arguments and create out directories
2. Loops over your list of genes and for each gene executes the script quality_check_reads_in_DADA2.R, which:
3. Grabs all the files in your data directory that match the patterns -p (default: _R1.fastq) and -q (default: _R2.fastq). This means it will grab all your data files but not anything else you might have in there.
4. Uses the dada2 function "plotQualityProfile" to plot quality graphs for -k (default: 24) samples. The samples are randomly pulled using the base-R `sample` function. You can use all your samples, but runtimes are long and the program doesn't do a great job of splitting the graphs up across multiple pages (potentially, a better functionality for higher k values is coming soon).
5. Writes those plots (one set for R1s, one set for R2s) to a pdf in a folder called "reports" inside your project directory.

The end result is a pair of PDFs for each gene called your_project_geneN_R[1-2]s_quality_profiles.pdf. You can then look at these and make decisions about your filtering parameters, I recommend following the dada2 tutorial to better understand the filtering options. Use these PDFs to fill out your params_file for step 4.

### Step Four: Filter and BLAST your reads

**Scripts**
step_4_wrapper.sh
- step_4_by_pident_local.sh
- step_4_start.sh
- step_4_by_score_local.sh
- step_4_by_pident_remote.sh
- filter_dada_part1.R
- filter_dada_part2.R
- step_3_p1_make_filt_parameters.sh

**What These Scripts Do**
Step 4 does the following, in order:
1. Parses arguments and creates out directories
2. Parses the parameter file and builds a filtering script in R with the appropriate arguments and settings.
3. Executes the newly-written filter_and_process_dada.R script to filter the data, which does the following:
  4. Grabs your data files based on your -p and -q patterns
  5. Filters the data using dada2's `FilterAndTrim` function, using your params_file for settings.
  6. Subsamples 50 randomly-selected samples down to 100k reads. This speeds up the error rate estimation process as it was not really built for huge #s of reads/sample, and is very robust at <100k reads.
  7. Performs the dada2 function  `LearnErrors` for the R1s and R2s to learn error rates.
  8. Performs the dada2 function `dada`, [the central algorithm of the dada2 pipeline](https://www.nature.com/articles/nmeth.3869#methods). This reduces samples down to unique ASV (amplicon sequence variants) within the learned error limits.
  9. Permforms the dada2 function `mergePairs` which merges the R1 and R2 unique reads into a single sequence
  10. Performs the dada2 function `removeBimeraDenovo` with the parameters `method="consensus"` to remove chimeric sequences.
  11. Produces an outfile in the "reports" directory tracking how many reads were retained for each sample at each of steps 5-10. This will only include samples that had >0 reads retained after filtering.
  12. Produces a raw ASV table with read counts for each sample at each ASV and ouptuts it in your results_tables directory in your project folder.
  13. Filters ASVs per sample: sets a relative read abundance (RRA) cutoff, default is 0.001 (0.1%). For each sample, multiplies total reads in that sample by the RRA cutoff to get the real_cutoff (so for a sample with 10,000 reads total across any number of ASVs, the real_cutoff would be 10,000*0.001=10). Outputs a file for each sample called sample_seqs.txt in your `gene_dada_out` directory in your project folder. The file contains all sequences with greater than real_cutoff reads. It is a two-column comma-separated file, column 1 contains sequences, column 2 contains reads.

The end result of this part of each of the scripts are a raw ASV table (rows are all valid, non-chimeric ASVs, no filtering, columns are samples, values are read counts), and a per-sample sample_seqs.txt file that contains all ASVs

The program then creates a combined_ASV.fasta file, the input for your BLAST search. It does this by combining every sample's seqs.txt file and extracting the unique sequences. It assigns a unique header name (seq_00001 through seq_99999) to each ASV.

What happens next in the scripts depends on your inputs. There are a few options that define the way the program makes decisions:

**What the Arguments Do**

*cutoff=97: This is an important cutoff. This is the percent identity (pident or %ident in BLAST) that the program will use as its cutoff for a "good" hit. This cutoff should depend on your system and your gene, but is generally between 97-99. This value **must* be an integer.

return_low=TRUE/FALSE: If not hits above cutoff% pident are available, should the program return the top hit? If TRUE, returns best hit according to the "score_pident" argument. If FALSE, returns "No Hit" for that sequence.

*local=TRUE/FALSE: If you select TRUE, you use your local reference database to perform BLAST and assess hits. This is the recommended method. It is faster and more tailored to your system. You will likely obtain higher resolution and more accurate assignments with a well-curated local reference database than with remote BLAST. If you select FALSE, you use the `-remote` function of `ncbi-blast` to perform remote BLAST rather than using a local database.

score_pident="bitscore"/"pident": Depending on which option you pick, the program will evaluate BLAST hits either by assessing the bitscore or the percent identity (https://ravilabio.info/notes/bioinformatics/e-value-bitscore.html). If "bitscore", then the "top BLAST hits" are considered to be the hit(s) with the highest bitscore (ie, equally high bitscores). If "pident", then the top BLAST hits are the hits with the top percent identity, regardless of length. In both cases the program then checks the highest scoring hit(s) are against the cutoff and if `return_low` is set to FALSE, hits below that cutoff are discarded.

**What the programs do next**
  14. Reads your user input and runs either local or remote BLAST on your combined ASV fasta. It returns a tab-formatted out-results file that includes the length, percent identity, and bitscore of every hit for every ASV. 
  15. creates a list of no-hits, as BLAST does not return data for sequences with no hit. 
  16. Adds the species and taxid to each BLAST hit.
  17. Makes a list of your sequences from your input fasta, and loops over that list. For each sequence, it:
    18. Pulls every matching hit from the BLAST output
    19. Assesses it to determine the highest scoring hit(s) (see the previous section for a description of how the program "scores" the hits)
    20. Checks how many species are present in the highest scoring hit(s). If it is only one, the program stores the species and taxonomic information for that sequence as the "best hit". If the top BLAST hit(s) include multiple species, the program then fetches the taxonomy tree for each species among the highest scoring hits and walks up the tree, checking at each level (genus, family, class, order, phylum) whether there is taxonomic agreement among the best hits. If no agreement is found at the phylum level, the sequence is recorded as a no-hit. Otherwise, the lowest commonly-shared taxon is reported as the "best hit" for that sequence.
  21. Once top BLAST hits have been selected for all ASVs, the program formats a best_hits.txt output file with one row per ASV.
  22. Formats and writes the taxatable, a tab-separated txt file (parsable by excel) that contains a row for every sample/ASV that reports the sample, the ASV, the reads, and the identity, score, and taxonomic information of the best hit for that ASV.

# FAQ
**What if I want my database to include every organism available for X gene?**

This pipeline is for limited-taxa reference databases or remote querying of NCBI. It is not designed to manage the curation of a large database that includes, say, all inverts with COI sequences in NCBI. If you want a larger database, consider using another reference database building tool, for instance [RESCRIPt](https://github.com/bokulich-lab/RESCRIPt) or [bcDatabaser](https://bcdatabaser.molecular.eco/). This pipeline is being actively developed to better improve performance with other databases. Any FASTA file that contains unique sequences each with a valid NCBI accession will work in step 2 of the pipeline, but runtimes for step 2 will be very long for huge databases (>1000 taxa).


**What if I have a few of my own in-house sequences I want to add to my database?**

That's fine, but it will add some tedioius work upfront. First, check if the species you are including have taxIDS assigned in NCBI. All taxa with at least one sequence in the NCBI nucleotide or sra database have a taxid, and you can look it up here: [here](https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi). Put your homebrewed sequences in a file called extras_gene1_sequences.fasta with the header format ">unique_identifier species name taxid=unique_taxid" and put it in your reference_database directory. If your species has a taxID assigned in NCBI, use that as the unique_taxid. Otherwise, you can assign one yourself, I recommend a series of letters and numbers separated by underscores. Do not use a strictly numeric value, you will almost definitely accidentally assign a real taxID to your species. You can assign any unique_identifier to the sequence (to replace the accession number), but it must not include spaces or special characters and each sequence must be uniquely named. Next, also in your reference database directory, install and run ncbitax2lin. Then, unzip and modify the ncbitax2lin file. Add the unique_taxids you assigned to the end of the file and add taxonomic information for each species following the comma-separated format of the ncbitax2lin file (described in the header of the file). This is tedious to do by hand, so we really recommend submitting your homebrew sequences to NCBI if they're long enough!

**What if I want to use remote BLAST to query all of NCBI?**

Sure, that's fine, start the pipeline at step 3 and then run the step_4_remote_blast_option.sh or step_4_remote_wrapper.sh (for head nodes and SLURM clusters, respectively). Expect the BLAST search to take several hours depending on the number of sequences, and make sure you have a good internet connection. If running on the standalone servers/head nodes of RCDS, I recommend using `screen` for all steps, but especially step 4, so the run will continue even if your pipe breaks.

**Should I use the "bitscore" or "pident" method of assessing BLAST hits?**
See the section above for the explanation about what bitscore and % identity are, and how the program uses them. The answer overall depends on your locus, its resolution, and the taxa in your ref database. For most datasets, bitscore is probably the best option. However, because the program only considers the hits with the single highest bitscore value, many nearly-equally-good hits may be hiding behind slightly-different bitscores. These bitscores may be different because of query cover differences of only a couple of base pairs. For high-quality datasets and loci with generally high resolution, choosing the pident option is generally going to be more conservative and recover fewer species-level hits than choosing bitscore. I encourage you to try both options and compare the results - overall they will likely be very similar, the differences in pident and bitscore will be seen mostly in edge cases.

**Why doesn't this pipeline use e-value instead of bitscore or pident?**
E-values help asses likelihood of homology and are affected by database size. The conserved regions of metabarcoding primers almost guarantee homology above any standard e-value cutoff, and these databases are generally small, so e-values are high across the board.

**I want to do 16S bacterial metagenomics. Is this pipeline right for me?** Probably not. There is an immense, incredibly well-maintained array of resources for 16S microbial amplicon sequencing, some commercial and some open source, that will be infinitely better than this pipeline.

#Installation and Usage
## Inputs and Preparation

### Inputs
This pipeline requires the following inputs:
- Demultiplexed, appropriately-named paired-end metabarcoding data in fastq (or fastq.gz) format, with adapters trimmed. Each gene or primer set should be analyzed separately, and should be in a different folder called gene1, gene2...geneN.
- A parameters file (see template) with your filtering parameters. You don't need this until step 3, so you can fill it out after you've seen the quality report on your reads. It must be in exactly the format as the template.
- The output of a single run of ncbitax2lin. This script will attempt to install and run ncbitax2lin if it can't find a taxonomy file. If you are having issues with permissions and install, you can run it yourself by copying the commands from [this readme](https://github.com/zyxue/ncbitax2lin), the script assumes a name of "ncbi_lineages_[date_of_utcnow].csv.gz" but is agnostic to the date, and you can provide a name if preferred.
- A taxa file containing the scientific names of your target taxa, one per line, with the header "taxname", if you want to run step 1 and query NCBI for your fasta.
- A file with lists of common terms for your target genes. One gene per column, as many permutations on the gene as you'd like, one per line (ie, "Cytochrome Oxidase I", "COI", "COX1"). See the example files for a template. The columns should have a header (which you can repeat in the body) that is a short, human-readable name of the gene that contains no spaces, slashes, quote marks or other special characters. For example, instead of heading your column "Cytochrome Oxidase I", head it "COI". This file can serve as your "genelist" file for steps 2 through 4 as well.

Data is often demultiplexed by the sequencing service company at no (or minor) cost. However, if your data has not been demultiplexed, we recommend using either [fastq-multx](https://github.com/brwnj/fastq-multx) or the demux-by-name function of [BBMap](https://github.com/BioInfoTools/BBMap). Some demultiplexers (fastq-multx, for instance) do automatic trimming, others do not, which is why this pipeline only accepts pre-trimmed data, although you have the option of using dada2's `trim` function to trim N bases from left/right. We recommend using [cutadapt](https://cutadapt.readthedocs.io/en/stable/) to trim. 

### Software and Installation

You will need
- Command-line git (pre-installed on UI RCDS Cluster, UCD FARM and Barbera)
- Command-Line BLAST [Installation Instructions Here](https://www.ncbi.nlm.nih.gov/books/NBK279690/) (pre-installed on UI RCDS Cluster, UCD FARM and Barbera)
- NCBItax2lin: [see here](https://github.com/zyxue/ncbitax2lin) used to add taxonomic information to the local reference database. Requires pip to install.
- gnu-sed (gsed) [Installation Instructions Here](https://formulae.brew.sh/formula/gnu-sed). Nearly all unix-based computing clusters use this as the default (including UI RCDS, FARM and Barbera), but it needs to be installed on a mac and aliased to 'sed'. It can be easily installed with [homebrew](https://brew.sh/)
- R: This pipeline was optimized for R version 4.2.3 
  Packages required:
  - optparse
  - tidyverse
  - lubridate
  - rentrez
  - stringr
  - bioconductor
  - dada2 (see dadad2 installation instructions [here](https://benjjneb.github.io/dada2/dada-installation.html) for bioconductor installation of dada2, may need to remove the "version=" flag to install.
**Providing an NCBI API Key** 
We strongly recommend providing Step 1 with an NCBI API key if you want to run the program on the SLURM nodes (recommended). Getting an NCBI account is free: [sign up here](https://account.ncbi.nlm.nih.gov/signup/). Once logged in, click your name in the upper right hand corner of the screen, click "Account Settings" and scroll down to the button that says "generate API key". Copy this key into the step_1_wrapper.sh script and you'll be good to go.


**If you are using the UI RCDS cluster:**
There is a good tutorial for installing r packages here: [link to RCDS](https://www.hpc.uidaho.edu/compute/Applications/R.html)
However, you need to make one change to the tutorial: you must load the newest R module (`module load R/4.2.3`). 
You can install your R packages wherever you'd like, but must supply the path to the scripts in the wrapper.

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

If you simply provide a project directory (folder) containing the step 1 input files and also containing a directory called gene1 (through ...geneN) that contains your fastas, the pipeline will create the rest of the directories for you within that project directory.

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
**Usage**
Edit the step_1_wrapper.sh file using a text editor like `nano` or `emacs` (NOT MS WORD) to include your relevant project names and paths. To specify your current working directory  as your project directory, you can use $PWD. For example, in the wrapper, `dirr=$PWD` would tell the script to use your current directory as the project directory. Similarly, `genelist=$PWD/genelist` will tell the script to look for a file called "genelist" in your current directory.

Options in the wrapper:
dirr=$PWD :full path to your_project directory. If you're running this code from inside this folder (recommended), you can use $PWD here.
db_dirr=reference_database :name (not path) to the out directory for your reference database. Default is "reference_database". Path is ${dirr}/reference_database.
genelist=$PWD/genelist :full path to your tab-delimited list of gene/primer set search terms.
taxlist=$PWD/taxlist :full path to your list of taxa. Space-delimited, single-column, header must be called "taxnames"
prefix="your_project" :name for your outfiles. Outfiles will generally be named "your_project_gene1" for each gene/primer set.
retmax=20 :maximum number of NCBI sequences to return. Recommended value is <100. 
rlib="~/Rpackages" :your R package library path
key="your_ncbi_api_key" : If you are running this on a cluster, we strongly recommend getting an NCBI API key because rate-limiting doesn't always work. They are free. Simply sign up (https://account.ncbi.nlm.nih.gov/signup) and navigate to >Account Settings. Click "generate API key" and put that code into this space.

Then run
`sbatch step_1_wrapper.sh` (for the SLURM node fortyfive)

or

`bash step_1_wrapper.sh` (for the head nodes like zaphod)


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

**Usage**
Fill in the blanks in the step_2_wrapper.sh script, including the path to your project directory and the **name** (not path) of the reference_database directory you chose. This directory must be nested inside your project directory, and it must contain at least 1 fasta file with the name scheme *_gene1_sequences.fasta for gene1...geneN.

dirr=/path/to/project/directory (can be $PWD if you're running the script from inside the project directory)
db_dirr=reference_database :name, not path, of database directory]
genelist=/path/to/genelist :the path to the file containing the gene terms
prefix="your_project" :a name for the output files, ex. "your_project", should be the same name as step 1.

Then run:

`sbatch step_2_wrapper.sh` (SLURM nodes)

or

`bash step_2_wrapper.sh` (zaphod, marvin, petunia, whale head nodes)

  
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
genelist=$PWD/genelist :path to your genelist, a tab-separated single-line list of genes/primer sets that correspond to your data folders.
prefix="your_project" :for ease of use, this should match the name you gave in Step 1.
pattern1="_R1.fastq" :the pattern for the end of your R1 files. Default is "_R1.fastq".
pattern2="_R2.fastq" :the pattern for the ened of your R2 files. Default is "_R2.fastq".
num_graphs=24 :change this based on how many quality profiles you want to look at. Adding more adds computational time to this script.`
rlib="~/Rpackages" :path to your R packages

then run:

`sbatch step_3_wrapper.sh` 

or 

`bash step_3_wrapper.sh`

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

dirr=$PWD :path to the your_project directory. If you are running this code from inside that directory (recommended), you can put $PWD here.
prefix=your_project :must use the same name as step_3
db_dir=reference_database :name, not path, of database directory (path will be ${dir}/${db_dir}). Default is 'reference_database'
pattern1="_R1.fastq" :default is _R1.fastq, leave blank for default
pattern2="_R2.fastq" :default is _R2.fastq, leave blank for default
genelist=$PWD/genelist :path to genelist
params_file="params_file" :name of params_file prefix. All parameter filenames should start with this prefix and end with _gene1,_gene2...geneN for each gene/primer set in the genelist. Default is params_file
localdat= :only put something here if your local database is not named "yourproject_gene_reference", which is the default for steps 1 and 2.
rlib="~/Rpackages" :path to your Rpackages
cutoff=98 :an integer cutoff (inclusive) for BLAST hits. hits below this are ignored if return_low is set to FALSE. 
return_low=TRUE Options: TRUE or FALSE. Set to TRUE to return the best available BLAST hit regardless of percent identity. Set to FALSE to discard hits below your cutoff. 
score="pident" Options: "pident" or "bitscore" only.
local="local" Options: "local" or "remote" only. Pick "local" to use your ref database to perform local BLAST or "remote" for remote BLAST (currently only supports the "pident" method).

then run:

`sbatch step_4_wrapper.sh` (fortyfive SLURM node)

or

`bash step_4_wrapper.sh` (zaphod, marvin, whale, petunia head nodes)


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

