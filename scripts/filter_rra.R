
args = commandArgs(trailingOnly=TRUE)
library(tidyverse,lib=args[1])
setwd(args[2])
full <- read.delim(args[3])
cutoff=as.numeric(args[5])
filt <- full %>% group_by(sample) %>% filter(reads >= (sum(reads)*cutoff))
outname=paste0(args[4],"_filtered_seqs.txt")
write_delim(filt, file=outname, quote="none", col_names=FALSE, delim = "\t")
