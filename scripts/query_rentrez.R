#!/usr/bin/R/
##2022 Shannon Kieran Blair##
##a script to query genbank for barcodes matching a list of taxa##

args = commandArgs(trailingOnly=TRUE)

library(rentrez, lib=args[2])
library(lubridate, lib=args[2])
library(tidyverse, lib=args[2])
library(stringr, lib=args[2])
library(optparse, lib=args[2])
#import your taxa and barcodes##
##INPUT FILE: tab-delimited list of genes, with header. Genes in columns. Rows are synonymous search terms. Taxlist: single column of taxnames, optional single column (tab delimited) for taxid

option_list <- list(make_option(c('-t','--taxfile'), action='store', type='character', default="example_taxfile.txt", help='path to taxfile, a single-column, space separated list of scientific names with the column header "taxname"'),
                    make_option(c('-g','--gene'), action='store', type='character', default="example_gene_list.txt", help='path to your tab-separated list of genes, one gene per column, one name per row (ie, COI, COX1'),
                    make_option(c('-d','--directory'), action='store', type='character', default="~", help='full path to output directory'),
                    make_option(c('-n', '--name'), action='store', type='character', default='your_project', help='prefix for naming outfiles'),
                    make_option(c('-z','--rlib'), action='store', type='character', default="~/Rpackages", help='r library path'),
                    make_option(c('-k','--key'), action='store', type='character', default="416858084f0b7a195302d9a4cff9627e5308", help='your NCBI api key'),
                    make_option(c('-r', '--retmax'), action='store', type='integer', default=10, help='how many matching fastas to return. Recommend <100')
)
opt <- parse_args(OptionParser(option_list = option_list))
set_entrez_key(opt$key)
Sys.getenv("ENTREZ_KEY")

##taxa format: taxid and taxname. Can use either one but need to update colnames to reflect data##
setwd(opt$directory)
prefix <- opt$name ##pick a name for your project. This name will be stuck on your filenames.
taxlist <-  as.data.frame(read_delim(opt$taxfile, delim="\t",show_col_types = FALSE))
genelist <- as.data.frame(read_delim(opt$gene,delim="\t",na="",show_col_types=FALSE))
genes <- colnames(genelist)
#set the loop parameters#
n <- ncol(genelist)
z <- nrow(taxlist)

##this loop first loops over each barcode, then loops over each taxa. 
##So each taxa is searched at each barcode. Stores hit data and the first FASTA result from each search. 
##Probaboly a much faster way to functionalize and apply this. 500 taxa and 3 barcodes takes about 50 minutes on my machine.
for(y in 1:n){
gene_terms <- str_c(na.omit(genelist[,y]), collapse=" OR ") ##I picked some common alternate names/spellings for each barcode. So the code loops over each column, which represents a barcode, and then creates a search term that includes all the alternate spellings I could think of.##
i=1
for (i in 1:z) {
  write(paste0("doing taxa ",i," out of ",z), stdout())
  Sys.sleep(0.1) 
  results <- entrez_search(db="nucleotide", term=paste0(taxlist$taxname[i],"[ORGN] AND (",gene_terms,")"),retmax=10) ##search the taxid and the barcode##
  if (length(results$ids) > 0) { ##if there are any results##
    avail_seq <- "yes"
    num_avail <- results$count ##number of matching results##
    Sys.sleep(0.1)
    summs <- entrez_summary(db="nucleotide", id=results$ids) ##summary of first result (it's so slow to grab all, plus it's tough to do anything with that much data)
    taxid <- unique(extract_from_esummary(summs, "taxid"))[1]
    accession <- extract_from_esummary(summs, "caption")
    Sys.sleep(0.1)
    temp.fasta <- entrez_fetch(db="nucleotide", id=results$ids, rettype="fasta")
    #temp.fasta <- sub("\n",paste0(" taxid=",taxid,"\n"),temp.fasta) #this adds the taxid to the first sequence. Useful later if you don't want to use taxid lookup to add id to seqs. But taxid lookup way is better and this causes weird behavior so for now i have commented it out#
    temp.data <- data.frame("taxname" = taxlist$taxname[i], "taxid" = as.character(taxid), "avail_seq"=avail_seq, "num_records"=num_avail) ##build individual dataframe.
    assign(paste0(taxlist$taxname[i],"_",colnames(genelist[y])),temp.data)
    genus=str_split_fixed(taxlist$taxname[i]," ",2)[1]
    species=str_split_fixed(taxlist$taxname[i]," ",2)[2]
    write(temp.fasta,file=paste0(prefix,"_",genus,"_",species,"_",colnames(genelist[y]),"_sequences.fasta"))
    rm(temp.data)
    rm(temp.fasta)
    rm(summs)
    } else {
    newr <- entrez_search(db="taxonomy", term=taxlist$taxname[i], retmax=1)
    if (length(newr$ids > 0)) {
      newsumms <- entrez_summary(db="taxonomy", id=newr$ids)
      taxid <- as.character(newsumms[7])
    } else {
      taxid <- "NA"
    }
    avail_seq <- "no"
    accession <- "NA"
    fasta_seq <- "NA"
    temp.data <- data.frame("taxname" = taxlist$taxname[i], "taxid" = taxid, "avail_seq"=avail_seq, "num_records"=0)
    assign(paste0(taxlist$taxname[i],"_",colnames(genelist[y])),temp.data)
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
if (length(genes) = 1 ) { 
write.table(paste0(colnames(genelist[y]),"_all"),file=paste0(prefix,"_database_taxa_summary.txt"), row.names=FALSE, col.names=TRUE, sep="\t",quote=FALSE,eol="\n")
} else {
list_of_all2 <-  mget(ls(pattern ="_all"))
##your whole list: each taxID has a row for each gene that includes the accession # and the number of records##
dfAby <- list_of_all2 %>% reduce(inner_join, by='taxname',na_matches = "na")
all_tax_out <- data.frame("taxname"=dfAby$taxname, "taxid"=dfAby$taxid.x)
all_tax_all_genes <- bind_rows(list_of_all2)
all_tax_avail_wide <- pivot_wider(all_tax_all_genes, id_cols="taxname", names_from="gene", values_from="avail_seq", names_prefix="avail_seq_")
all_tax_num_wide <- pivot_wider(all_tax_all_genes, id_cols="taxname", names_from="gene", values_from="num_records", names_prefix="num_records_")
##another table you might want, just a list of which taxa have seq for which barcodes##
all_tax1 <- inner_join(all_tax_out,all_tax_avail_wide,by="taxname")
rm(all_tax_out)
all_tax_out <- inner_join(all_tax1, all_tax_num_wide,by="taxname")
as <- "avail_seq_"
cols1 <- paste0(as,unlist(genes))
all_tax_out <- unite(all_tax_out, "avail_temp", all_of(cols1), sep=" ", remove=FALSE)
all_tax_out$avail_any <- ifelse(grepl("yes", all_tax_out$avail_temp), "yes", "no")
all_tax_out <- all_tax_out[,-3]
rm(all_tax1)
rm(list_of_all2)
rm(all_tax_avail_wide)
rm(all_tax_all_genes)
write.table(all_tax_out, file=paste0(prefix,"_database_taxa_summary.txt"), row.names=FALSE, col.names=TRUE, sep="\t",quote=FALSE,eol="\n")
}
