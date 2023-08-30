args=commandArgs(trailingOnly=TRUE)
library(lubridate, lib=args[1])
library(tidyverse, lib=args[1])

setwd(args[6])
df1=read.delim(args[2], header=FALSE)
colnames(df1) <- c("sequence","seqnum","local_identity","local_species","local_taxid","local_phylum","local_class","local_order","local_family","local_genus","local_bitscore","local_num_spec_in_best_hit","local_all_spec_in_best_hit")
df2=read.delim(args[3], header=FALSE)
colnames(df2) <- c("seqnum","remote_identity","remote_species","remote_taxid","remote_phylum","remote_class","remote_order","remote_family","remote_genus","remote_bitscore","remote_num_spec_in_best_hit","remote_all_spec_in_best_hit")

df3 <- inner_join(df1,df2,by="seqnum")

fname=paste0(args[4],"_",args[5],"best_blast_hits.txt")
write_delim(df3, fname, delim="\t",quote="none")
