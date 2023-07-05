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
		    make_option(c('-k','--num_graphs'), action='store', type='integer', default=24, help='prefix for naming your outfiles. default is 24') 
            )
opt <- parse_args(OptionParser(option_list = option_list))

setwd(opt$directory) ##set a directory that contains your filenames
path <- paste0(opt$directory,"/",opt$gene)  ##this is the path, likely it will be identical to your working directory.
#list.files(path) ##make a list of all the files in your path
write("plotting R1 quality profiles, this may take some time.",stdout())
fnFs <- sort(list.files(path, pattern =opt$pattern1)) ##CHANGE THIS to match your F/R or R1/R2 naming structure.
fnRs <- sort(list.files(path, pattern =opt$pattern2)) ##see above

setwd(path)
sample.names <- sapply(strsplit(basename(fnFs),"_"),`[`,1) ##change this pattern too - it will currently cut the basename off
samp_fnFs <- sample(fnFs, size=opt$num_graphs)
samp_fnRs <- sample(fnRs, size=opt$num_graphs)
pdf(file=paste0(opt$name,"_",opt$gene,"_R1s_quality_profiles.pdf"), onefile=TRUE)
print(plotQualityProfile(samp_fnFs)) #I look at a lot of these, but only about 30-50 at a time. So the default should be about 1:50, and then I run it again at 51:100, etc.
dev.off()
write("plotting R2 quality profiles, this may take some time.",stdout())
pdf(file=paste0(opt$name,"_",opt$gene,"_R2s_quality_profiles.pdf"), onefile=TRUE)
print(plotQualityProfile(samp_fnRs)) #I look at a lot of these, but only about 30-50 at a time. So the default should be about 1:50, and then I run it again at 51:100, etc.
dev.off()

