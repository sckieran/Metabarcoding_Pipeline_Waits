# Metabarcoding_Pipeline_Waits
A rough pipeline and tutorial to analyze metabarcoding data using dada2. This explains my system. This project is theoretically able to run on the UIdaho RCDS cluster and on a mac running Monterey. However, it is intended to be an analysis guide for students working on their own projects, rather than a button-press way to receive your output. This is a teaching tool, so it is pretty wordy. Please read through it to understand all the decision points in this analysis.

To use this pipeline on the cluster, you will likely need to alter the Rscripts to include the line .libPaths("~/path/to/your/R/library") after the shebang line, and that path will need to direct the script to a folder containing your installed R packages. To install packages onto the UI RCDS cluster, follow the tutorial [here:](https://www.hpc.uidaho.edu/compute/Applications/R.html)

This pipeline assumes reasonably good resolution of locus data and is designed as a "first pass" at your metabarcoding data, not a finished product. **This pipeline grabs the top BLAST hit for each unique sequence. It does not assess multiple BLAST hits.** There is an optional extension script that will attempt to resolve identically-good BLAST hits. This extension script was built for a specific project where the resolution of the locus was well-understood and there were relatively few possible target taxa. It is extremely naive and basic. If your resolution is low (ie, you have many disparate taxa with perfect or near-perfect matches at your locus), or if you don't understand the resolution of your locus, consider using a Bayesian taxonomic assignment algorithm like the ones found in DADA2 or QIIME instead.

**This pipeline is not designed for microbial metagenomics at 16S**. I say this because there is an immense, incredibly well-maintained array of resources for 16S microbial amplicon sequencing, some commercial and some open source, that will be infinitely better than this pipeline.

## Inputs and Installation

### Inputs
This pipeline requires the following inputs:
- Demultiplexed, appropriately-named paired-end metabarcoding data in fastq (or fastq.gz) format, with adapters trimmed. Each gene or primer set should be analyzed separately, and should be in different folders.
- A tab-separated parameters file containing: a list of genes/loci/primer sets and the expected length of each amplicon.
- The output of a single run of ncbitax2lin. You can copy paste the commands from [this readme](https://github.com/zyxue/ncbitax2lin), the script assumes a name of "ncbi_lineages_[date_of_utcnow].csv.gz" but is agnostic to the date, and you can provide a name if preferred.

Additionally, to run the local database building tool, you need:
- A taxa file containing the scientific names of your target taxa, one per line, with the header "taxname"
- A file with lists of common terms for your target genes. One gene per column, as many permutations on the gene as you'd like, one per line (ie, "Cytochrome Oxidase I", "COI", "COX1"). See the example files for a template. The columns should have a header (which you can repeat in the body) that is ta short, human-readable name of the gene that contains no spaces, slashes, quote marks or other special characters. For example, instead of heading your column "Cytochrome Oxidase I", head it "COI".

Data is often demultiplexed by the sequencing service company at no (or minor) cost. However, if your data has not been demultiplexed, we recommend using either [fastq-multx](https://github.com/brwnj/fastq-multx) or the demux-by-name function of [BBMap](https://github.com/BioInfoTools/BBMap). Some demultiplexers (fastq-multx, for instance) do automatic trimming, others do not. The dual-indexed fusion primers often used in metabarcoding may require extra trimming of the overhangs. We recommend using [cutadapt](https://cutadapt.readthedocs.io/en/stable/) to trim. 

### Software and Installation

The pipeline is relatively light on software. R package management will probably be the most intensive thing.

- Python 3 is needed to run many programs. It is usually installed when needed (ie, cutadapt), but having a python install on your Mac is a good thing, especially now that there is no Python 2 installed natively on new Macs.
- Command-Line BLAST [Installation Instructions Here](https://www.ncbi.nlm.nih.gov/books/NBK279690/)
- NCBItax2lin: [see here](https://github.com/zyxue/ncbitax2lin) used to add taxonomic information to the local reference database. Requires pip to install.
- gnu-sed (gsed) [Installation Instructions Here](https://formulae.brew.sh/formula/gnu-sed). Most computing clusters use this as the default, but it needs to be installed on a mac. It can be easily installed with [homebrew](https://brew.sh/)
- R: This pipeline was optimized for R version 4.1.2 -- Bird Hippie.
- Packages required:
  - optparse
  - tidyverse
  - rentrez
  - dada2
  - stringr

### Before you Begin
The pipeline assumes all internal R scripts are in the same folder as the shell scripts. 

## Optional Step 0: Build a Local Reference Database

We recommend a local reference database to improve the accuracy of detections. Future versions of this pipeline may also include a method to format this library so it can be used with DADA2's taxonomic assignment algorithm. This script uses [rentrez](https://cran.r-project.org/web/packages/rentrez/vignettes/rentrez_tutorial.html), an excellent R package for querying NCBI. However, you can skip this step and use remote BLAST instead. Building a reference library may take several hours, but the actual searches are much faster than remote BLAST.

## Things to Consider
- This script assumes a relatively small reference library, <500 taxa at <5 genes. This script is rate limited to avoid overloading the NCBI API, but can still run into odd errors when internet connections are unstable, and doesn't always play nicely with VPNs. If you're getting weird errors like "unexpected EOF" or "http 400", try disconnecting any VPNs or checking firewalls. I am generally unable to help you troubleshoot these issues.
- We recommend that you include common contaminants in your database for better error detection. At minimum, we recommend *Homo sapiens* for any genes that might amplify vertebrates and *Cyprinus carpio* for all projects. As of mid-2023, the *Cyprinus carpio* genome is full of Illumina adapters and therefore if you have a lot of primer dimer in your reads, or if you fail to trim your adapters, you will have many hits to this species. I assume it will get a new genome at some point, so this may not be true forever. 
- This script can search at any taxonomic level. For species, include the binomial name. For subspecies, check NCBI for the correct name of the subspecies you're interested in, or use taxids. For genera and above, use the single-word taxon name.
- There is an extension to the R script that can check for sequences within a genus for any species that has no available relevant seq in genbank, to add congeners and therefore accuracy to the database. By default, it is commented out. To use it, simply remove the comment character (#) inside the query_rentrez.R script, starting at the section labeled "search no-hits by genus" and going to the end of the file.
- This script doesn't check for the wrong organism in the search results. This is because in most cases, multiple taxa returned from a search are subspecies, or because the search term was for a genus or higher taxonomic rank. However, there is a version of this script called step_0_make_local_database_by_taxid.sh that you can use instead if you are worried about sister taxa slipping into your search results. The input for that script is a taxfile containing a list of taxids, which you are responsible for supplying. In theory, you can use [this NCBI tool](https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi) to look up taxids. **The taxid version of this script only accepts species/subspecies as input, not genera or any higher taxonomic rank, and will check to ensure that only the target taxid is included in the output files.**
- 
### Inputs

- A taxa file containing the scientific names of your target taxa, one per line, with the header "taxname"
- A tab-delimited file with lists of genes with "gene terms". Each column corresponds to a gene, each row is different ways to describe that gene. This just provides the best search radius. So 12S would include "12S" "12S Ribosomal RNA" "12S Mitochondrial Sequence". The columns can be of different lengths. We recommend 2-3 terms per gene for the most common metabarcoding loci (12S, COI, ITS, 16S)


### Usage

**Arguments/Options**
The script requires four command line arguments:
* -n a name for the output files, ex. "your_project"
* -t the path to the file containing a list of taxa
* -g the path to the file containing the gene terms
* -d the path to the output directory

**Optional Arguments:**
* -r (retmax) maximum number of search records to return. Default 10. Setting this value very high (>50) may cause problems with the rentrez search and will probably add off-target sequences (wrong gene more likely than wrong organism) to your database. I recommend getting an NCBI API key for large databases, see [the rentrez tutorial](https://cran.r-project.org/web/packages/rentrez/vignettes/rentrez_tutorial.html#rate-limiting-and-api-keys) for more info.
* -c combine: choice of "comb" or "sep" (default). If "comb" is invoked, the fastas from all genes will be combined into a single database. This is not recommended. Sequences from different genes/loci must be processed separately (because they will almost certainly be different lengths)
  
bash step_0_make_local_database.sh -n your_name -t /path/to/taxfile -g /path/to/genefile -d /path/to/outputdir -r 10 -c sep

## Output

In your output directory:
- A tab-delimited summary file with information about each taxon including the taxid, number of available sequences on NCBI (at any locus), whether sequences are available for any of your target loci, and a yes/no for each locus, called "your_project_database_taxa_summary.txt"
- A list of taxa that had no hits for your relevant genes in NCBI, called "your_project_database_taxa_no_hits.txt"
- A fasta file for each gene/locus called your_project_gene_database_seqs.fasta
- Directories for each gene containing the fasta files for each taxa separately, which can be deleted if you're done with them. 
- A searchable local NCBI reference database for each gene. These consists of 10 files, each with the name "your_project_gene_reference.[suffix]".


##End of current content##


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


