##2022 Shannon Kieran Blair##
##a script to query genbank for barcodes matching a list of taxa##

setwd("/your/data/directory")  ##change this to your directory##
library(rentrez)
library(tidyverse)
library(tidyr)
library(stringr)
library(purrr)
#import your taxa and barcodes##
##INPUT FILE: tab-delimited list of genes, with header. Genes in columns. Rows are synonymous search terms. Taxlis: single column of taxnames, optional single column (tab delimited) for taxid

##taxa format: taxid and taxname. Can use either one but need to update colnames to reflect data##
taxlist <- read.delim("list_of_taxa.txt", header=TRUE)
genelist <- read.delim("list_of_genes.txt",na.strings="")

n <- ncol(genelist)
z <- nrow(taxlist)
##this loop first loops over each barcode, then loops over each taxa. 
##So each taxa is searched at each barcode. Stores hit data and the first FASTA result from each search. 
##Probaboly a much faster way to functionalize and apply this. 500 taxa and 3 barcodes takes about 50 minutes on my machine.
for(y in 1:n){
gene_terms <- str_c(na.omit(genelist[,y]), collapse=" OR ") ##I picked some common alternate names/spellings for each barcode. So the code loops over each column, which represents a barcode, and then creates a search term that includes all the alternate spellings I could think of.##
i=1
for (i in 1:z) {
  results <- entrez_search(db="nucleotide", term=paste0(taxlist$taxname[i],"[ORGN] AND (",gene_terms,")"),retmax=10) ##search the taxid and the barcode##
  #print(results)
  if (length(results$ids) > 0) { ##if there are any results##
    avail_seq <- "yes"
    num_avail <- results$count ##number of matching results##
    summs <- entrez_summary(db="nucleotide", id=results$ids[1]) ##summary of first result (it's so slow to grab all, plus it's tough to do anything with that much data)
    taxid <- summs[[10]]
    accession <- summs[[3]]
    temp.fasta <- entrez_fetch(db="nucleotide", id=results$ids, rettype="fasta")
    temp.data <- data.frame("taxname" = taxlist$taxname[i], "taxid" = as.character(taxid), "avail_seq"=avail_seq, "num_records"=num_avail,"accession"=accession) ##build individual dataframe.
    assign(paste0(taxlist$taxname[i],"_",colnames(genelist[y])),temp.data)
    write(temp.fasta,file=paste0("neil_local_",taxlist$taxname[i],"_",colnames(genelist[y]),"sequences.fasta"))
    #rm(results)
    rm(temp.data)
    rm(temp.fasta)
    rm(summs)
  } else {
  avail_seq <- "no"
  accession <- "NA"
  fasta_seq <- "NA"
  taxid <- "NA"
  temp.data <- data.frame("taxname" = taxlist$taxname[i], "taxid" = taxid, "avail_seq"=avail_seq, "num_records"=0,"accession"=accession, "fasta"=fasta_seq)
  assign(paste0(taxlist$taxname[i],"_",colnames(genelist[y])),temp.data)
  #rm(results)
  rm(temp.data)
}
}
##aggregate into a single table and remove all those extra data frames##
list_of_gene <-  mget(ls(pattern = paste0("_",colnames(genelist[y]))))
temp_gene_all <- bind_rows(list_of_gene)
temp_gene_all$gene <- colnames(genelist[y]) ##add a column with barcode ID
assign(paste0(colnames(genelist[y]),"_all"),temp_gene_all)
rm(list = ls(pattern = paste0("_",colnames(genelist[y]),"$"))) ##clean up data
rm(temp_gene_all)
}

list_of_all2 <-  mget(ls(pattern ="_all"))
##your whole list: each taxID has a row for each gene (so should be 3x3) that includes the accession # and FASTA for the first hit for each one, and the number of records##
all_tax_all_genes <- bind_rows(list_of_all2)
write.table(all_tax_all_genes, file="your_project_localref_full_taxa_list.txt", row.names = FALSE,col.names = TRUE,sep="\t",quote = FALSE, eol="\n")
##probably the table you really want, just a list of which taxa have seq for which barcodes##
all_tax_wide <- pivot_wider(all_tax_all_genes, id_cols="taxname", names_from="gene", values_from="avail_seq")
write.table(all_tax_wide, file="your_project_localref_summary.txt", row.names=FALSE, col.names=TRUE, sep="\t",quote=FALSE,eol="\n")
all_pos_hits <- all_tax_all_genes %>% filter(avail_seq=="yes")
seqlist <- data.frame("taxname"=all_pos_hits$taxname, "seq"=all_pos_hits$fasta)
write.table(seqlist, file="your_project_localref_sequences.txt", row.names=FALSE, col.names=FALSE, quote=FALSE, eol="\n")

##exploring our no-hits##

#get list of searches with no hits##
##this is not a generalized solution. It is specific to having 3 barcodes. Writing a more generalized solution is on the to-do list. It's not complicated, but I got lazy trying to figure out the fastest count method##
no_hits <- all_tax_all_genes %>% filter(avail_seq == "no")
##let's filter out what's missing all three barcodes vs. what is only missing one or two barcodes##
dups_test <- no_hits[which(duplicated(no_hits$taxname)),] ##this is a list of those missing at least 2 barcodes##
no_hits_any_barcode <- dups_test[which(duplicated(dups_test$taxname)),] ##this is a list of those missing all 3 barcodes
new_taxlist <- dups_test$taxname

##looping again. Yes, probably a way to functionize this and apply it. This grabs the search results for any sequences associated with that organism (and we switch to name over taxid just in case), then it prints how many sequences are available. If it's mostly 0 or 1, that's a good sign our searches are working and maybe can highlight a few worth checking out by hand.#
for (q in 1:length(new_taxlist)) {
  ser <- entrez_search(db="nucleotide", term=paste0(new_taxlist[q],"[ORGN]"))
  summ <- entrez_summary(db="nucleotide", id=ser$ids)
  assign(paste0(new_taxlist[q],"_fullrecord"),summ)
  total_seqs <- length(ser$ids)
  temp.data <- data.frame("taxname" = new_taxlist[q], "total_seqs" =total_seqs)
  assign(paste0(new_taxlist[q],"_seqnum"),temp.data)
  rm(temp.data)
  
}

#rm(list = ls(pattern = "_fullrecord$"))
##aggregate list of taxa with no barcodes and how many seqs they have total in genbank, for further examination later##
list_of_seqnums <- mget(ls(pattern ="_seqnum"))
no_hit_total_seqs <- bind_rows(list_of_seqnums)
write.table(no_hit_total_seqs,file="your_project_localref_no_hits_seq_results.txt", row.names=FALSE, col.names=TRUE, sep="\t",quote=FALSE,eol="\n")
rm(list = ls(pattern = "_seqnum"))




