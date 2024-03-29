args = commandArgs(trailingOnly=TRUE)
library(lubridate, lib=args[7])
library(tidyverse, lib=args[7])

setwd(args[2])

ts <- read.delim(args[1])
un <- unique(ts$sample)

for (i in 1:length(un)) {
t1 <- ts %>% filter(sample == un[i] & identity > args[6])
test <- ts %>% filter(sample == un[i] & identity > args[6]) %>% group_by(taxa) %>% summarise(sum(reads))
cutoff <- sum(test$`sum(reads)`) * as.numeric(args[5])
passing <- test %>% filter(`sum(reads)` > cutoff)
t1 <- t1[t1$taxa %in% passing$taxa,]
assign(paste0(i,"_t1"),t1)
}
rm(t1)
list <- mget(ls(pattern = "_t1"))
new_ts <- bind_rows(list)
rm(list=ls(pattern = "_t1"))
fil=paste0(args[2],"/",args[3],"_",args[4],"_filtered_taxatable.txt")
write_delim(new_ts,file=fil,quote="none",delim="\t")
