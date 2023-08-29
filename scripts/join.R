args = 
library(lubridate, lib=args[1])
library(tidyverse, lib=args[1])

setwd(args[6])
df1=read.delim(args[2], header=FALSE)
df2=read.delim(args[3], header=FALSE)

df3 <- inner_join(df1,df2,by=V2)

fname=paste0(args[4],"_",args[5],"best_blast_hits.txt")
write_delim(df3, fname, delim="\t",quote="none")
