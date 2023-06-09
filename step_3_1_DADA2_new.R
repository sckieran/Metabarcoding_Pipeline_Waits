library(dada2) ##load your liibrary
library(tidyverse)
setwd("/path/to/your/raw/files") ##set a directory that contains your filenames
path <- "/path/to/your/raw/files"  ##this is the path, likely it will be identical to your working directory.
list.files(path) ##make a list of all the files in your path

fnFs <- sort(list.files(path, pattern ="_R1.fastq")) ##CHANGE THIS to match your F/R or R1/R2 naming structure.
fnRs <- sort(list.files(path, pattern ="_R2.fastq")) ##see above

sample.names <- sapply(strsplit(basename(fnFs),"_"),`[`,1) ##change this pattern too - it will currently cut the basename off
plotQualityProfile(fnFs[1:10]) #I look at a lot of these, but only about 30-50 at a time. So the default should be about 1:50, and then I run it again at 51:100, etc.
plotQualityProfile(fnRs[1:25]) #if you want to look really closely you can also just look at the first 2 or 4.
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq")) #set filtered output names
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names


###this filters and truncates your reads. Set your trunclength based on your overlap size ###
###and where the quality in your reads drop. See the DADA2 tutorial for more guidance on filtering parameters.###

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,truncLen=c(200,200), ##maxN must be 0 but all other filters can be changed
                     maxN=0, maxEE=c(2,2), rm.phix=TRUE,
                     compress=FALSE, multithread=TRUE) # On Windows set multithread=FALSE 
head(out)

###
##These steps learn the error rates, plot them for a gut check, and then preps for the ASV calling step###
###

errF <- learnErrors(filtFs, multithread=TRUE) #learning error rates
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)


##these are the actual ASV calling steps, they call ASVs and then merge F/R###

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
head(mergers[[1]])

###this produces a table of merged ASVs and looks at how big it is##
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

##this step removes chimeras, ie falsely merged ASVs##

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

### this process, which currently doesn't work if you have samples where 0 reads pass filters###
### shows you how many reads were discarded at each step. There shouldn't be any step with a big "cliff" of read fallout.###
###check the dada2 tutorial for troubleshooting ideas if you have big cliffs of fallout (ie, if it happens at mergers then you over-truncated etc.)###

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.csv(track, "austin_p1thru6_its_sampletracking.txt", row.names = FALSE, quote = FALSE)

###these steps format and write your ASV tables to be used with BLAST or other tools###

seqt <- as.data.frame(t(seqtab.nochim))  #transpose your ASV file so rows are seqs and columns are samples
x <- ncol(seqt)
#setwd("~/Documents/Documents_Local/student_metabarcoding/austen/dada_output/trnl_p56")
##this loop takes each sample, filters the ASVs that occur at a rate of >0.01% of total reads (for that sample) and writes ASV and read count to a file##
for (i in 1:x) {
  no <- sum(seqt[,i])
  cutoff <- 0.001 * no ##filtering threshold, change this to change threshold##
  col1 <- as.data.frame(seqt) %>% filter(seqt[,i] > cutoff)
  temp_col <- col1[,i]
  temp_col <- as.data.frame(temp_col)
  rownames(temp_col) <- rownames(col1)
  write.csv(x=temp_col,file=paste0(colnames(seqt[i]),"_seqs.txt"), row.names=TRUE, col.names=FALSE, quote=FALSE)
}

##this turns the sequence table into a dataframe, filters the data to only ASVs with read counts >0.01% of the entire dataset##
##then it writes the whole dataset to a file to be turned into a fasta for BLAST##
##NOTE: THIS FILTERING IS NOT EQUAL TO THE FILTERING ABOVE##
##THIS WILL FILTER MORE READS BECAUSE IT FILTERS ON % OF OVERALL READS RATHER THAN PER-SAMPLE##

seqtd <- as.data.frame(seqt)
nos <-  colSums(seqtd)
no <- sum(nos)
cutoff <- 0.00001 * no ##filtering threshold, change this to change threshold##
seqtd$sums <- rowSums(seqtd)
seqt_filt <- seqtd %>% filter(sums > cutoff)
write.csv(seqt_filt, "austin_allsamples_trnl_p5p6_nofilter_seqlist.txt", row.names=TRUE, quote=FALSE)

###other filtering options##
##filter on <1% present in ANY sample (much looser, even looser than above), by using max rather than rowsums, check dplyr cheatsheet for details##
##write a few lines to count reads in negative controls (or just look), use highest val in negs as cutoff##
##look for contamination in the positive control, use that read count as cutoff##
##take some average between neg and pos control##
##use flat cutoff, like 100 reads##
##literature is not settled on this, so add a few lines to look at what you filtered:
seqout_filtered <- seqtd %>% filter(sums < cutoff)
write.csv(seqout_filtered, "austin_allsamples_ITS_filtered_at_001_excludedseqs.txt", row.names=TRUE, quote=FALSE)

#filter what makes sense#
