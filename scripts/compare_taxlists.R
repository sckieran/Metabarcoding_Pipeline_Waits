
args = commandArgs(trailingOnly=TRUE)
library(tidyverse, lib=args[1])
setwd(args[2])
local <- read.delim(args[3])
remote <- read.delim(args[4])

remote_t <- data.frame("sequence"=remote$sequence,"seqnum"=remote$seqnum,"remote_pident"=remote$identity,"remote_taxa"=remote$species)
local_t <- data.frame("sequence"=local$sequence, "local_pident"=local$identity, "local_taxa"=local$species)
 
local_filt <- local_t[local_t$sequence %in% remote$sequence,]
comp_hits <- inner_join(local_filt,remote_t, by="sequence")
write.table(comp_hits, file=paste0(args[5],"_",args[6],"_comparative_BLAST_hits.txt"), sep="\t", quote=FALSE, row.names=FALSE)
