##Visualize FALCON results from top hits file
library(ggplot2)
library(taxonomizr)
args=commandArgs(trailingOnly=TRUE)
targ=args[[1]];
taxdmp = args[[2]];


top = read.delim(targ, sep="\t", header=FALSE)
#parse accession numbers, ref length, and percent similarity to separate table

r = 1;
accmatrix= matrix(ncol=ncol(top)+1, nrow=nrow(top)*5);
for (i in 1:nrow(top)){
  split = strsplit(as.character(top[i,4]), "_")
  grsplit = grep("[.]", split[[1]]);
  for(n in grsplit){
    accmatrix[r,1:7]= unlist(c((split[[1]][n]), top[i,])); #Warning fasta header returned here as factor ID rather than character string
    r = r+1;
  }
}
accmatrix=na.omit(accmatrix)



taxaId<-accessionToTaxa(accmatrix[,1],paste(taxdmp, "accessionTaxa.sql", sep='/'))
taxaNodes<-read.nodes(paste(taxdmp, 'nodes.dmp', sep = '/')) 
taxaNames<-read.names(paste(taxdmp, 'names.dmp', sep = '/')) 

taxtable= getTaxonomy(taxaId,taxaNodes,taxaNames)

taxtable = cbind(taxtable, accmatrix)
taxtable = data.frame(taxtable, row.names=NULL)
colnames(taxtable) = c('superkingdom', 
                       'phylum', 
                       'class', 
                       'order', 
                       'family', 
                       'genus', 
                       'species', 
                       'accNo', 
                       'hitNo', 
                       'length', 
                       'similarity', 
                       'nameFactor', 
                       'dbStart', 
                       'dbEnd')

targpath = strsplit(targ, "/")
targbase = strsplit(targpath[[1]][length(targpath[[1]])], "[.]")
 targbase = targbase[[1]][length(targbase[[1]])-1];
# dir.create(targbase)
# setwd(targbase)
phylum = stats::na.omit(plyr::count(taxtable$phylum))

kingdom = stats::na.omit(plyr::count(taxtable$superkingdom))
png(filename = "FALCONkingdom.png", height = 7, width = 11, units = "in", res = 600)
kingdom <- kingdom[order(-kingdom$freq),]
kingdom$x <- factor(kingdom$x, levels = kingdom$x[order(-kingdom$freq)])
kingdom$freq = kingdom$freq/sum(kingdom$freq)*100
ggplot(data = kingdom, aes_string(x = colnames(kingdom)[1], y = "freq")) + geom_col() +  labs(x = "SuperKingdom") + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + ggtitle(targbase)
dev.off()
 
family = stats::na.omit(plyr::count(taxtable$family))
png(filename = "FALCONfamily.png", height = 7, width = 11, units = "in", res = 600)
family <- family[order(-family$freq),]
family$x <- factor(family$x, levels = family$x[order(-family$freq)])
family$freq = family$freq/sum(family$freq)*100
ggplot(data = family[1:50,], aes_string(x = colnames(family)[1], y = "freq")) + geom_col() +  labs(x = "Family") + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + ggtitle(targbase)
dev.off()
# 
# fam2 = aggregate(as.numeric(as.character(taxtable$length)) ~ taxtable$family, FUN='sum')
# colnames(fam2) = c("family", "length")
# png(filename = "FALCONfamilylen.png", height = 7, width = 11, units = "in", res = 600)
# fam2 <- fam2[order(-fam2$length),]
# fam2$family <- factor(fam2$family, levels = fam2$family[order(-fam2$length)])
# ggplot(data = fam2[1:50,], aes_string(x = colnames(fam2)[1], y = "length")) + geom_col() +  labs(x = "Family") + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + ggtitle(targbase)
# dev.off()

genus = stats::na.omit(plyr::count(taxtable$genus))
png(filename = "FALCONgenus.png", height = 7, width = 11, units = "in", res = 600)
genus <- genus[order(-genus$freq),]
genus$x <- factor(genus$x, levels = genus$x[order(-genus$freq)])
genus$freq = genus$freq/sum(genus$freq)*100
ggplot(data = genus[1:50,], aes_string(x = colnames(genus)[1], y = "freq")) + geom_col() + labs(x = "Genus") + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + ggtitle(targbase)
dev.off()

#plants only
pfamily = family[which(family$x %in% subset(taxtable, taxtable$phylum=='Streptophyta')$family),]
png(filename = "FALCONpfamily.png", height = 7, width = 11, units = "in", res = 600)
ggplot(data = pfamily[1:50,], aes_string(x = colnames(pfamily)[1], y = "freq")) + geom_col() + labs(x = "Plant genera") + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + ggtitle(targbase)
dev.off()

pgenus = genus[which(genus$x %in% subset(taxtable, taxtable$phylum=='Streptophyta')$genus),]
png(filename = "FALCONpgenus.png", height = 7, width = 11, units = "in", res = 600)
ggplot(data = pgenus[1:50,], aes_string(x = colnames(pgenus)[1], y = "freq")) + geom_col() + labs(x = "Plant genera") + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + ggtitle(targbase)
dev.off()

write.table(taxtable, file = 'FTax.Results', sep = '\t')
write.table(kingdom, file = 'Fsuperkingdom', sep = '\t')
write.table(family, file = 'Fsuperfamily', sep = '\t')
write.table(phylum, file = 'Fphylum', sep = '\t')
write.table(order, file = 'Forder', sep = '\t')
write.table(genus, file = 'Fgenus', sep = '\t')







