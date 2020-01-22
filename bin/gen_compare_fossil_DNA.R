# read data for macrofossils and aDNA
root = "~/NeotomaSeq"
files = list.files(root, recursive=T, full.names = T, pattern="*.ResultsAfter")

# Macrofossils
macro = read.table("./COR_macrofossils.csv", sep=',', header =T, stringsAsFactors = F)

colnames(macro) = macro[1,]
macro=macro[-1,]
head(macro)

#compile taxa by site ##NOTE: Taxa are taken from the 'names' column where we keep only the first word. This is either Genus or Family
#add family collectiono
n=1
taxa = list()
fams = list()
sites = vector()
for(i in 5:ncol(macro)){
  names = vector()
  for(z in macro$Name[which(macro[,i]>0)]){
    names = append(names, unique(strsplit(z, ' ')[[1]][1]))
  }
  famvec = vector()
  for(z in macro$Family[which(macro[,i]>0)]) {
    famvec=append(famvec, unique(z));
  }
  taxa[[n]] = names
  fams[[n]] = gsub('Amaranthaceae', replacement = 'Chenopodiaceae', x = famvec)
  sites[n] = colnames(macro)[i]
  n=n+1
}
print(taxa)

## Read the aDNA genus files
dnagenfiles = files
dnagenera = list()
dnasites = vector()
dnagencount = list()

for(i in dnagenfiles){
  print(i)
  if(grepl('Pr', i)){next}
  if(grepl('GC100', i)){next}
  if(grepl('TS211A', i)){next}
  nom = strsplit(i, '/Taxa.ResultsAfter')[[1]][1]
  nom2 = strsplit(nom, '[/]')[[1]]
  name = nom2[length(nom2)]
  dnagen = data.table::fread(i, sep = '\t')
  #get plant hits ONLY
  dnasub = dnagen[dnagen$phylum == 'Streptophyta',]
  sumreads = sum(dnagen$numUniqueReads)
  aggr = aggregate(dnasub$numUniqueReads, by=list(dnasub$genus), FUN=sum )
  subs = aggr[aggr$x>(0.01*sumreads),]
  dnasites[name] = name
  dnagenera[[name]] = as.character(subs$Group.1)
  dnagencount[[name]] = subs$x
}

## Read the aDNA family files
dnafamfiles = files
dnafamily = list()
dnafamsites = vector()
dnafamcount = list()
for(i in dnafamfiles){
  print(i)
  if(grepl('Pr', i)){next}
  if(grepl('GC100', i)){next}
  if(grepl('TS211A', i)){next}
  nom = strsplit(i, '/Taxa.ResultsAfter')[[1]][1]
  nom2 = strsplit(nom, '[/]')[[1]]
  name = nom2[length(nom2)]
  dnafam = data.table::fread(i, sep = '\t')
  #get plant hits ONLY
  dnasub = dnafam[dnafam$phylum == 'Streptophyta',]
  #sumreads = sum(dnasub$numUniqueReads)
  aggr = aggregate(dnasub$numUniqueReads, by=list(dnasub$family), FUN=sum )
  sumreads = sum(aggr$x)
  subs = aggr[aggr$x>(0.01*sumreads),]
  dnafamsites[name] = name
  dnafamily[[name]] = as.character(subs$Group.1)
  dnafamcount[[name]] = subs$x
}

##FAMILY PLOTTING
allfoss=unique(c(unlist(fams), unlist(dnafamily)))
nr = length(allfoss)
coll = matrix(nrow=nr, ncol = 21)
xx = 1
columns = vector()
for(nn in 1:length(dnasites)){
  q=dnasites[[nn]]
  
  presite = strsplit(q, '[-]')[[1]][1]
  r = grep(presite, sites)
  if(length(r)==0){next}
  columns[[xx]] = paste(sites[[r]], '.foss', sep='')
  
  coll[which(allfoss %in% fams[[r]]),xx] = 1
  
  xx = xx + 1
  
  columns[[xx]] = paste(dnasites[[nn]], '.dnagen', sep ='')
  coll[which(allfoss %in% dnafamily[[nn]]),xx] = dnafamcount[[nn]]/max(dnafamcount[[nn]])
  xx= xx + 1
}

coll= as.data.frame(coll)
colnames(coll) = columns
print(coll)
coll = cbind(allfoss, coll)
#plot using pheatmap or ggplot2 tile
library(ggplot2)
library(plyr)
library(scales)
library(tidyr)
library(dplyr)
coll.m <- melt(coll)

coll.m = coll.m %>% 
  separate(variable, c('site'), sep = '[-.]', remove = F, extra='drop') %>% 
  filter(!grepl("Pr", variable)) %>% filter(!is.na(site))

(macrofam <- ggplot(coll.m, aes(variable, allfoss)) + 
    geom_tile(aes(fill = value), colour = "black") +
    scale_fill_gradient2(low = "white", mid = muted("blue"),
                         high = "navyblue", midpoint = 0.4, space = "Lab",
                         na.value = "grey95", guide = "colourbar", aesthetics = "fill") +
    #scale_fill_gradient(low = "purple", high = "navyblue", na.value = 'white') +
    facet_grid(~site, scales = 'free_x', space='free') +
    theme_linedraw() + 
    xlab('') + ylab('') +
    labs(fill='Rel. Abundance') +
    theme(
      axis.text = element_text(size = 6),
      axis.text.x  = element_text(angle = 45, size =8, hjust=1),
      axis.text.y  = element_text(angle = 45, size =8, hjust=1),
      axis.title = element_text(size = 6),
      legend.text = element_text(size = 6),
      legend.position="top"
    ))


png(filename = 'macrofam.png', height = 7.5, width = 5.75, units = 'in', res=400)
macrofam
dev.off()





##GENUS PLOTTING
allfoss=unique(c(unlist(taxa), unlist(dnagenera)))
allfoss = allfoss[grep('aceae', allfoss, invert=TRUE)]
nr = length(allfoss)
coll = matrix(nrow=nr, ncol = 21)
xx = 1
columns = vector()
for(nn in 1:length(dnasites)){
  q=dnasites[[nn]]
  
  presite = strsplit(q, '[-]')[[1]][1]
  r = grep(presite, sites)
  if(length(r)==0){next}
  columns[[xx]] = paste(sites[[r]], '.foss', sep='')
  
  coll[which(allfoss %in% taxa[[r]]),xx] = 1
  
  xx = xx + 1
  
  columns[[xx]] = paste(dnasites[[nn]], '.dnagen', sep ='')
  coll[which(allfoss %in% dnagenera[[nn]]),xx] = dnagencount[[nn]]/max(dnagencount[[nn]])
  xx= xx + 1
}

coll= as.data.frame(coll)
colnames(coll) = columns
print(coll)
coll = cbind(allfoss, coll)
#plot using ggplot2 tile
library(ggplot2)
library(plyr)
library(scales)
library(tidyr)
library(dplyr)
coll.m <- melt(coll)

coll.m = coll.m %>% 
  separate(variable, c('site'), sep = '[-.]', remove = F, extra='drop') %>% 
  filter(!grepl("Pr", variable)) %>% filter(!is.na(site))

(macrogen <- ggplot(coll.m, aes(variable, allfoss)) + 
    geom_tile(aes(fill = value), colour = "black") +
    scale_fill_gradient2(low = "white", mid = muted("blue"),
                         high = "navyblue", midpoint = 0.4, space = "Lab",
                         na.value = "grey95", guide = "colourbar", aesthetics = "fill") +
    #scale_fill_gradient(low = "purple", high = "navyblue", na.value = 'white') +
    facet_grid(~site, scales = 'free_x', space='free') +
    xlab('') + ylab('') +
    labs(fill='Rel. Abundance') +
    theme_linedraw() + 
    theme(
      axis.text = element_text(size = 6),
      axis.text.x  = element_text(angle = 45, size =8, hjust=1),
      axis.text.y  = element_text(angle = 45, size =8, hjust=1),
      axis.title = element_text(size = 6),
      legend.text = element_text(size = 6),
      legend.position="top"
    ))


png(filename = 'macrogen.png', height = 7.5, width = 5.75, units = 'in', res=400)
macrogen
dev.off()


library(ggpubr)
gsa = ggarrange(macrofam, macrogen, ncol=2, nrow=1, labels="AUTO")
ggsave(filename='macrocomp.pdf', plot=gsa, device=NULL, width = 9.25, heigh=6.5, dpi=600)
ggsave(filename='macrocomp.png', plot=gsa, device=NULL, width = 9.25, heigh=6.5, dpi=500)


