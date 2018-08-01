#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

library(stringr)
library(taxonomizr)
taxdmp = args[1] #taxdump folder
all_hits = read.delim2(file = args[2], header = FALSE, sep = "\t") #filtered BLast hits
#Pull out accession numbers
subject <- all_hits[which(all_hits$V2=="Subject"),]
ID = strsplit(as.character(subject$V3), ":")
ID[[1]][1]
ACCID = vector()
for(i in 1:length(ID)){
  ACCID[i] = ID[[i]][1]
}
ACCID = unique(ACCID)
#Get taxonomic classifications from accession numbers
taxaId<-accessionToTaxa(ACCID,paste(taxdmp, "accessionTaxa.sql", sep='/'))
taxaNodes<-read.nodes(paste(taxdmp, 'nodes.dmp', sep = '/')) 
taxaNames<-read.names(paste(taxdmp, 'names.dmp', sep = '/')) 

taxtable= getTaxonomy(taxaId,taxaNodes,taxaNames)
