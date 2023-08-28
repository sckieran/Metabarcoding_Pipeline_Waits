
args = commandArgs(trailingOnly=TRUE)
library(lubridate,lib=args[5])
library(tidyverse,lib=args[5])
setwd(args[1])
full <- read.delim(args[2])
cutoff=as.numeric(args[4])
colnames(full) <-c("seq","reads")
filt <- full %>% filter(reads >= (sum(reads)*cutoff) | reads >= 500)
outname=paste0(args[3],"_filtered_seqs.txt")
write_delim(filt, file=outname, quote="none", col_names=FALSE, delim = "\t")

