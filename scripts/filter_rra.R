
args = commandArgs(trailingOnly=TRUE)
library(lubridate, lib=args[1])
library(tidyverse,lib=args[1])
setwd(args[2])
full <- read.delim(args[3])
cutoff=as.numeric(args[5])
colnames(full) <-c("seq","reads")
filt <- full %>% filter(reads >= (sum(reads)*cutoff) | reads >= 500)
outname=paste0(args[4],"_filtered_seqs.txt")
write_delim(filt, file=outname, quote="none", col_names=FALSE, delim = "\t")

