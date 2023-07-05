args = commandArgs(trailingOnly=TRUE)
library(dada2, lib=args[2]) ##load your library
library(lubridate, lib=args[2])
library(tidyverse, lib=args[2])
library(optparse, lib=args[2])

option_list <- list(make_option(c('-d','--directory'), action='store', type='character', default="~", help='path to the enclosing directory, which contains a folder called "gene" that contains your demultiplexed and trimmed fastqs'),
                    make_option(c('-g','--gene'), action='store', type='character', default="12S", help='path to your tab-separated list of genes, one gene per column, one name per row ie COI, COXI'),
                    make_option(c('-p','--pattern1'), action='store', type='character', default="_R1.fastq", help='pattern for naming your R1 files, default is _R1.fastq'),
                    make_option(c('-q','--pattern2'), action='store', type='character', default="_R2.fastq", help='pattern for naming your R2 files, default is _R2.fastq'),
                    make_option(c('-n','--name'), action='store', type='character', default="your_project", help='prefix for naming your outfiles. default is "your_project"'),
                    make_option(c('-e','--rracutoff'), action='store', type='numeric', default=0.0005, help='filter ASVs at readcounts < rracutoff.'),
		    make_option(c('-z','--rlib'), action='store', type='character', default="~/Rpackages", help='r library path'),
		    make_option(c('-m','--multithread'), action='store', type='logical', default="TRUE", help='TRUE or FALSE to multithread threadable functions')
)
opt <- parse_args(OptionParser(option_list = option_list))

setwd(opt$directory) ##set a directory that contains your filenames
path <- paste0(opt$directory,opt$gene)  ##this is the path, likely it will be identical to your working directory.

outpath <- paste0(opt$directory,opt$gene,"_dada_out")
subpath <- paste0(outpath,"/subsampled")
list.files(path) ##make a list of all the files in your path
write(paste0("path to data files is, ",path),stdout())
write(paste0("your naming pattern is ",opt$pattern1," for R1s and ",opt$pattern2," for R2s."),stdout())
fnFs <- sort(list.files(path, pattern=opt$pattern1)) ##CHANGE THIS to match your F/R or R1/R2 naming structure.
fnRs <- sort(list.files(path, pattern=opt$pattern2)) ##see above
setwd(path)
sample.names <- sapply(strsplit(basename(fnFs),opt$pattern1),`[`,1) ##change this pattern too - it will currently cut the basename off

filtFs <- file.path(outpath, "filtered", paste0(sample.names, "_F_filt.fastq")) #set filtered output names
filtRs <- file.path(outpath, "filtered", paste0(sample.names, "_R_filt.fastq"))
write(paste0("path to seq files is, ",outpath),stdout())
write(paste0("path to filtered data is, ",outpath,"/filtered"),stdout())
names(filtFs) <- sample.names
names(filtRs) <- sample.names

###this filters and truncates your reads. Set your trunclength based on your overlap size ###
###and where the quality in your reads drop. See the DADA2 tutorial for more guidance on filtering parameters.###
##to make things easier, after this step I go back and remove reads that have 0 samples passing filter, then re-run this script##

write("filtering data, this may take some time. You can watch progress in the filt directory you specified.", stdout())
