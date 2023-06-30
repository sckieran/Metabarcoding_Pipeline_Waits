
###
##These steps learn the error rates, plot them for a gut check, and then preps for the ASV calling step###
###

filtpath=paste0(outpath,"/filtered")  ##these steps re-call your filtered files. I did this because some files had 0 reads pass filter after filtration. This messes up your downstream QC though. I'm working on a better solution.
filtFs_list <- sort(list.files(filtpath, pattern ="_F_filt.fastq"))
filtRs_list <- sort(list.files(filtpath, pattern ="_R_filt.fastq"))
filtFs <- file.path(filtpath, filtFs_list)
filtRs <- file.path(filtpath, filtRs_list)

sample.names <- sapply(strsplit(basename(filtFs_list),"_F_filt.fastq"),`[`,1)
setwd(filtpath)
system("ls *_F_filt.fastq | shuf -n50 > sublist")
system("while read k; do base=$(echo ${k} | awk -F\"_F_filt.fastq\" '{print $1}'); head -n400000 ${k} > ../subsampled/${base}_sub_R1.fastq; head -n400000 ${base}_R_filt.fastq > ../subsampled/${base}_sub_R2.fastq; done < sublist")
subFs<- sort(list.files(subpath, pattern=paste0("_sub",opt$pattern1), full.names = TRUE))
subRs<- sort(list.files(subpath, pattern=paste0("_sub",opt$pattern2), full.names = TRUE))
write("learning error rates", stdout())

errF <- learnErrors(subFs, multithread=opt$multithread) #learning error rates
errR <- learnErrors(subRs, multithread=opt$multithread)


##performing the dada function###
write("now performing dada function")
dadaFs <- dada(filtFs, err=errF, multithread=opt$multithread)
dadaRs <- dada(filtRs, err=errR, multithread=opt$multithread)

#merging read pairs##
write("now merging read pairs", stdout())
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

###this produces a table of merged ASVs and looks at how big it is##
write("formatting your data")
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

##this step removes chimeras, ie falsely merged ASVs##
write("removing chimeras",stdout())
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=opt$multithread, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

### this process, which currently doesn't work if you have samples where 0 reads pass filters###
### shows you how many reads were discarded at each step. There shouldn't be any step with a big "cliff" of read fallout.###
###check the dada2 tutorial for troubleshooting ideas if you have big cliffs of fallout (ie, if it happens at mergers then you over-truncated etc.)###
write("now producing dada2 processing report")
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
write(track, file=paste0(opt$directory,"/reports/",opt$name,"_",opt$gene,"_dada_report.txt"),sep = "\t")

###these steps format and write your ASV tables to be used with BLAST or other tools###

write("filtering your samples at your designated relative read abundance value and producing sequence files")
seqt <- as.data.frame(t(seqtab.nochim))  #transpose your ASV file so rows are seqs and columns are samples
x <- ncol(seqt)
write.table(x=seqt, file=paste0(opt$directory,"/results_tables/",opt$name,"_",opt$gene,"_raw_ASV_table.txt"),quote=FALSE) 
##this loop takes each sample, filters the ASVs that occur at a rate of >0.01% of total reads (for that sample) and writes ASV and read count to a file##
for (i in 1:x) {
  no <- sum(seqt[,i])
  cutoff <- opt$rracutoff * no
  col1 <- as.data.frame(seqt) %>% filter(seqt[,i] > cutoff)
  temp_col <- col1[,i]
  temp_col <- as.data.frame(temp_col)
  rownames(temp_col) <- rownames(col1)
  write.table(x=temp_col,file=paste0(opt$directory,"/",opt$gene,"_dada_out/",colnames(seqt[i]),"_seqs.txt"),sep=",", row.names=TRUE, col.names=FALSE, quote=FALSE)
}
