#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
library(parallel)
library(CHNOSZ)
library(ggplot2)
library(dplyr)
library(taxonomizr)
#getNamesAndNodes()

nclus = 4
taxaNodes<-read.nodes.sql(paste(args[[2]], '/nodes.dmp', sep = '')) #takes time?
taxaNames<-read.names.sql(paste(args[[2]], '/names.dmp', sep = '')) 

#Table -> taxID and numreads
taxdump = args[2]

#Read the table into a cleaner file
cent.table = read.delim(args[[1]], header = FALSE, sep = "\t", skip = 1) #row 1 = bad headers

# chtax = function(x) {
# 
#   require(dplyr)
#   
#   cent.small = cent.table[x, ]
# 
#   #Create empty vectors
#   parentG=vector();
#   parentO=vector();
#   parentK=vector();
#   parentP=vector();
#   parentF=vector();
#   parentSK=vector();
# 
#   #Get all parents of cent.table
#   parentG = CHNOSZ::parent(cent.small$V2, taxdir= taxdump, rank="genus")
#   parentO = CHNOSZ::parent(cent.small$V2, taxdir= taxdump, rank="order")
#   parentK = CHNOSZ::parent(cent.small$V2, taxdir= taxdump, rank="kingdom")
#   parentP = CHNOSZ::parent(cent.small$V2, taxdir= taxdump, rank="phylum")
#   parentF = CHNOSZ::parent(cent.small$V2, taxdir= taxdump, rank="family")
#   parentSK = CHNOSZ::parent(cent.small$V2, taxdir= taxdump, rank="superkingdom")
# 
#   parentK[which(parentK == 1)] = parentSK[which(parentK == 1)]
#   
#   parentK[which(parentK == 2759)] = 1
#   
#   #megablast
#   parent.genera = CHNOSZ::sciname(parentG, taxdump)
#   parent.orders = CHNOSZ::sciname(parentO, taxdump)
#   parent.kingdoms = CHNOSZ::sciname(parentK, taxdump)
#   parent.phylums = CHNOSZ::sciname(parentP, taxdump)
#   parent.families = CHNOSZ::sciname(parentF, taxdump)
# 
#   #parent.genera and table
#   cent.small$V8<-NA
#   cent.small <- cent.small %>%
#     mutate(V8 = parent.genera)
# 
#   #parent.order and table
#   cent.small$V9<-NA
#   cent.small <- cent.small %>%
#     mutate(V9 = parent.orders)
#   
#   #parent.kingdom and table
#   cent.small$V10<-NA
#   cent.small <- cent.small %>%
#     mutate(V10 = parent.kingdoms)
#   
#   
#   #parent.phylum and table
#   cent.small$V11<-NA
#   cent.small <- cent.small %>%
#     mutate(V11 = parent.phylums)
#   
#   
#   #parent.family and table
#   cent.small$V12<-NA
#   cent.small <- cent.small %>%
#     mutate(V12 = parent.families)
#   
# 
#   #Order, genus, numreads, original reads table 
#   graph.table <- data.frame(cent.small$V10, cent.small$V11, cent.small$V9, cent.small$V12, cent.small$V8, cent.small$V1, cent.small$V2, cent.small$V3, cent.small$V4, cent.small$V5, cent.small$V6, cent.small$V7)
#   names(graph.table) <- c("Kingdoms", "Phylums", "Orders", "Families", "Genera", "Original_Name", "Tax_ID", "Tax_Rank", "Genome_Size", "Num_Reads", "Num_Unique_Reads", "Abundance")
#   
#   return(graph.table)
# }

chtax2 = function(x) {
  require(taxonomizr)
  cent.small = cent.table[x, ]
  taxtable= getTaxonomy(cent.small$V2,taxaNodes,taxaNames)
  
  return(taxtable)
}

do = seq(1, nrow(cent.table))
cl = makeCluster(nclus, type = "SOCK")
p = proc.time()
splits = clusterSplit(cl, do)
clusterExport(cl, list("cent.table", "taxdump", "taxaNodes", "taxaNames"))
par.table = parLapply(cl, splits, chtax2)

stopCluster(cl)

#Initialize final.table
final.table = par.table[[1]]

#Loop to create table formatted correctly
for(i in 2:length(par.table)){
  final.table = rbind(final.table, par.table[[i]])
}

write.table(final.table, file = "Taxa.ResultsBefore", quote = FALSE, sep = "\t", row.names = FALSE)

final.table <- cbind(final.table, cent.table)
names(final.table) <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species", "origName", "taxID", "taxRank", "genomeSize", "numReads", "numUniqueReads", "abundance")

#final.table <- filter(final.table, Kingdoms != "root", Phylums != "root", Orders != "root", Genera != "root", Families != "root" )

#Write table to a file
write.table(final.table, file = "Taxa.ResultsAfter", quote = FALSE, sep = "\t", row.names = FALSE)


#agregate (dplyr) -- which column to collapse to same value summaries numreads by sum function
summary.genera <- aggregate(apply(final.table[,c('numReads', 'numUniqueReads')], 2, as.numeric), by = list(final.table$genus), FUN = sum)  
summary.orders <- aggregate(apply(final.table[ ,c('numReads', 'numUniqueReads')], 2, as.numeric), by = list(final.table$order), FUN = sum)
summary.kingdoms <- aggregate(apply(final.table[ ,c('numReads', 'numUniqueReads')], 2, as.numeric), by = list(final.table$superkingdom), FUN = sum)
summary.phylums <- aggregate(apply(final.table[ ,c('numReads', 'numUniqueReads')], 2, as.numeric), by = list(final.table$phylum), FUN = sum)
summary.families <- aggregate(apply(final.table[,c('numReads', 'numUniqueReads')] , 2, as.numeric), by = list(final.table$family), FUN = sum)

#Sort and plot (bar chart)
png(filename = "GeneraPlot.png", height = 7, width = 11, units = "in", res = 600 )
summary.genera <- summary.genera[order(-summary.genera$numUniqueReads),]
summary.genera$Group.1 <- factor(summary.genera$Group.1, levels = summary.genera$Group.1[order(-summary.genera$numUniqueReads)])
ggplot(data = summary.genera[1:50,], aes_string(x = colnames(summary.genera)[1], y = 'numUniqueReads')) + geom_col() + scale_y_continuous(trans='log10') + labs(x = "Genera") + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
dev.off()

png(filename = "OrdersPlot.png", height = 7, width = 11, units = "in", res = 600)
summary.orders <- summary.orders[order(-summary.orders$numUniqueReads),]
summary.orders$Group.1 <- factor(summary.orders$Group.1, levels = summary.orders$Group.1[order(-summary.orders$numUniqueReads)])
ggplot(data = summary.orders[1:50,], aes_string(x = colnames(summary.orders)[1], y = 'numUniqueReads')) + geom_col() + scale_y_continuous(trans='log10') + labs(x = "Orders") + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
dev.off()

png(filename = "KingdomsPlot.png", height = 7, width = 11, units = "in", res = 600)
summary.kingdoms <- summary.kingdoms[order(-summary.kingdoms$numUniqueReads),]
summary.kingdoms$Group.1 <- factor(summary.kingdoms$Group.1, levels = summary.kingdoms$Group.1[order(-summary.kingdoms$numUniqueReads)])
ggplot(data = summary.kingdoms, aes_string(x = colnames(summary.kingdoms)[1], y = 'numUniqueReads')) + geom_col() + scale_y_continuous(trans='log10') + labs(x = "Kingdoms") + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
dev.off()

png(filename = "PhylumsPlot.png", height = 7, width = 11, units = "in", res = 600)
summary.phylums <- summary.phylums[order(-summary.phylums$numUniqueReads),]
summary.phylums$Group.1 <- factor(summary.phylums$Group.1, levels = summary.phylums$Group.1[order(-summary.phylums$numUniqueReads)])
ggplot(data = summary.phylums[1:50,], aes_string(x = colnames(summary.phylums)[1], y = 'numUniqueReads')) + geom_col() + scale_y_continuous(trans='log10') + labs(x = "Phylums") + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
dev.off()

png(filename = "FamiliesPlot.png", height = 7, width = 11, units = "in", res = 600)
summary.families <- summary.families[order(-summary.families$numUniqueReads),]
summary.families$Group.1 <- factor(summary.families$Group.1, levels = summary.families$Group.1[order(-summary.families$numUniqueReads)])
ggplot(data = summary.families[1:50,], aes_string(x = colnames(summary.families)[1], y = 'numUniqueReads')) + geom_col() + scale_y_continuous(trans='log10') + labs(x = "Families") + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
dev.off()

name <- strsplit(basename(args[[1]]), '[.]')

write.table(summary.genera, file =  paste(name[[1]][1], "Genus", sep = "."), quote = FALSE, sep = "\t", row.names = FALSE)
write.table(summary.orders, file = paste(name[[1]][1], "Order", sep = "."), quote = FALSE, sep = "\t", row.names = FALSE)
write.table(summary.kingdoms, file = paste(name[[1]][1], "Superkingdom", sep = "."), quote = FALSE, sep = "\t", row.names = FALSE)
write.table(summary.phylums, file = paste(name[[1]][1], "Phylum", sep = "."), quote = FALSE, sep = "\t", row.names = FALSE)
write.table(summary.families, file = paste(name[[1]][1], "Family", sep = "."), quote = FALSE, sep = "\t", row.names = FALSE)
