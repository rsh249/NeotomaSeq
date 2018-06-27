##Visualize FALCON results from top hits file
library(ggplot2)
library(taxonomizr)

setwd('~/FALCONVIS')
top = read.delim("top.csv", sep="\t", header=FALSE)
#parse accession numbers, ref length, and percent similarity to separate table

r = 1;
accmatrix= matrix(ncol=ncol(top)+1, nrow=nrow(top)*5);
for (i in 1:nrow(top)){
  split = strsplit(as.character(top[i,4]), "_")
  grsplit = grep("[.]", split[[1]]);
  for(n in grsplit){
    accmatrix[r,1:7]= unlist(c((split[[1]][n]), top[i,]))
    r = r+1;
  }
}
accmatrix=na.omit(accmatrix)


taxaId<-accessionToTaxa(accmatrix[,1],"accessionTaxa.sql")
taxaNodes<-read.nodes('nodes.dmp') 
taxaNames<-read.names('names.dmp')

taxtable= getTaxonomy(taxaId,taxaNodes,taxaNames)

taxtable = cbind(taxtable, accmatrix)
taxtable = data.frame(taxtable, row.names=NULL)

dir.create('./FALCONfigs')
setwd('FALCONfigs')


family = stats::na.omit(plyr::count(taxtable$family))
png(filename = "FALCONfamily.png", height = 7, width = 11, units = "in", res = 600)
family <- family[order(-family$freq),]
family$x <- factor(family$x, levels = family$x[order(-family$freq)])
ggplot(data = family[1:50,], aes_string(x = colnames(family)[1], y = "freq")) + geom_col() +  labs(x = "Family") + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
dev.off()

genus = stats::na.omit(plyr::count(taxtable$genus))
png(filename = "FALCONgenus.png", height = 7, width = 11, units = "in", res = 600)
genus <- genus[order(-genus$freq),]
genus$x <- factor(genus$x, levels = genus$x[order(-genus$freq)])
ggplot(data = genus[1:50,], aes_string(x = colnames(genus)[1], y = "freq")) + geom_col() + labs(x = "Genus") + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
dev.off()

#plants only
plants = subset(taxtable, taxtable$phylum=='Streptophyta')
pfamily = stats::na.omit(plyr::count(plants$family))
png(filename = "FALCONpfamily.png", height = 7, width = 11, units = "in", res = 600)
pfamily <- pfamily[order(-pfamily$freq),]
pfamily$x <- factor(pfamily$x, levels = pfamily$x[order(-pfamily$freq)])
ggplot(data = pfamily[1:50,], aes_string(x = colnames(pfamily)[1], y = "freq")) + geom_col() + labs(x = "Plant genera") + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
dev.off()

pgenus = stats::na.omit(plyr::count(plants$genus))
png(filename = "FALCONpgenus.png", height = 7, width = 11, units = "in", res = 600)
pgenus <- pgenus[order(-pgenus$freq),]
pgenus$x <- factor(pgenus$x, levels = pgenus$x[order(-pgenus$freq)])
ggplot(data = pgenus[1:50,], aes_string(x = colnames(pgenus)[1], y = "freq")) + geom_col() + labs(x = "Plant genera") + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
dev.off()



