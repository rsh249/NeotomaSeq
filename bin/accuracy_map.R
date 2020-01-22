# read data for:

# read data for aDNA
root = "~/NeotomaSeq"
files = list.files(root, recursive=T, full.names = T, pattern="*.ResultsAfter")



# Macrofossils

macro = read.table("./COR_macrofossils.csv", sep=',', header =T, stringsAsFactors = F)

colnames(macro) = macro[1,]
macro=macro[-1,]
head(macro)

#compile families by site
n=1
fams = list()
sites = vector()
for (i in 5:ncol(macro)) {
  fams[[n]] = (unique(macro$Family[which(macro[, i] > 0)]))
  sites[n] = colnames(macro)[i]
  n = n + 1
}
dna2macro_coll = list()
nn=1
for (cutoff in c(0.00001, 0.0001, 0.001, 0.01, 0.1)) {
  ## Read the aDNA family files
  dnafamfiles = files
  dnafamily = list()
  dnafamsites = vector()
  dnafamcount = list()
  for (i in dnafamfiles) {
    
    print(i)
    #if(grepl('Pr', i)){next}
    if(grepl('GC100', i)){next}
    if(grepl('TS211A', i)){next}
    nom = strsplit(i, '/Taxa.ResultsAfter')[[1]][1]
    nom2 = strsplit(nom, '[/]')[[1]]
    name = nom2[length(nom2)]
    dnafam = data.table::fread(i, sep = '\t')
    #get plant hits ONLY
    dnasub = dnafam[dnafam$phylum == 'Streptophyta', ]
    #sumreads = sum(dnasub$numUniqueReads)
    aggr = aggregate(dnasub$numUniqueReads,
                     by = list(dnasub$family),
                     FUN = sum)
    sumreads = sum(aggr$x)
    subs = aggr[aggr$x > (cutoff * sumreads), ]
    dnafamsites[name] = name
    dnafamily[[name]] = as.character(subs$Group.1)
    dnafamcount[[name]] = subs$x
  }
  
  
  #calc accuracy
  dna2macro = matrix(ncol = length(dnafamsites), nrow = 5)
  columns = vector()
  zz = 1
  
  for (i in 1:length(sites)) {
    b = grep(sites[i], dnafamsites)
    if (length(b) == 1) {
      #samplename = strsplit(dnasites[[b]], '/')[[1]][6]
      samplename = dnafamsites[[b]]
      columns[zz] = samplename
      summ = sum(fams[[i]] %in% dnafamily[[b]])
      dna2macro[1, zz] = summ
      dna2macro[2, zz] = summ / length(fams[[i]])
      dna2macro[3, zz] = summ / length(dnafamily[[b]])
      dna2macro[4, zz] = length(fams[[i]])
      dna2macro[5, zz] = length(dnafamily[[b]])
      zz = zz + 1
    } else {
      for (n in b) {
        samplename = dnafamsites[[n]]
        columns[zz] = samplename
        summ = sum(fams[[i]] %in% dnafamily[[n]])
        dna2macro[1, zz] = summ #total in common
        dna2macro[2, zz] = summ / length(fams[[i]]) #TPR_macro
        dna2macro[3, zz] = summ / length(dnafamily[[n]]) #TPR_DNA
        dna2macro[4, zz] = length(fams[[i]]) #totalmacro
        dna2macro[5, zz] = length(dnafamily[[n]]) #totalDNA
        zz = zz + 1
      }
    }
    
  }
  
  dna2macro = as.data.frame(dna2macro)
  colnames(dna2macro) = columns
  rownames(dna2macro) = c('intersect',
                          'TPR_macro',
                          'TPR_DNA',
                          'totalmacro',
                          'totalDNA')
  dna2macro = as.data.frame(t(dna2macro))
  dna2macro = cbind(columns, dna2macro)
  dna2macro_coll[[nn]] = cbind(dna2macro, rep(cutoff, nrow(dna2macro)))
  nn = nn+1
  p1 = ggplot(dna2macro) +
    geom_col(aes(x = columns, y = TPR_DNA)) +
    ylim(c(0, 1))
  p2 = ggplot(dna2macro) +
    geom_col(aes(x = columns, y = TPR_macro)) +
    ylim(c(0, 1))
  
  library(gridExtra)
  ggsave('TPR.png',
         grid.arrange(p1, p2))
  
  
  
}

accuracy_map = dna2macro_coll[[1]]
for(i in 2:length(dna2macro_coll)){
  accuracy_map = rbind(accuracy_map, dna2macro_coll[[i]])
}

colnames(accuracy_map)[7] = 'cutoff'
prec = ggplot() +
  geom_line(data=accuracy_map, aes(x=cutoff*100, y=(intersect)/(intersect + (totalDNA-intersect)))) +
  facet_grid(~columns) + scale_x_log10() + 
  theme_linedraw() + 
  ylab('Precision') + 
  xlab('') + 
  theme(
    axis.text = element_text(size = 8),
    axis.text.x  = element_text(angle = 45, size = 0, hjust =1),
    axis.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    strip.text = element_text(size=7)
  )


sens = ggplot() +
  geom_line(data=accuracy_map, aes(x=cutoff*100, y=(intersect)/(intersect + (totalmacro-intersect)))) +
  facet_grid(~columns) + scale_x_log10() + 
  theme_linedraw() + 
  ylab('Sensitivity') + 
  xlab('% of plant reads as threshold') + 
  theme(
    axis.text = element_text(size = 8),
    axis.text.x  = element_text(angle = 45, size = 8, hjust =1),
    axis.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    strip.text = element_text(size=7)
  )


FDR = ggplot() +
  geom_line(data=accuracy_map, aes(x=cutoff*100, y=(totalDNA - intersect)/((totalDNA - intersect) + intersect))) +
  facet_grid(~columns) + scale_x_log10() + 
  theme_linedraw() + 
  ylab('False Discover Rate') + 
  xlab('% of plant reads as threshold') + 
  theme(
    axis.text = element_text(size = 8),
    axis.text.x  = element_text(angle = 45, size = 8, hjust =1),
    axis.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    strip.text = element_text(size=7)
  )
FDR
options(scipen=999, digits=1) #shut off scientific notation

library(ggpubr)
gsa = ggarrange(prec, sens, ncol=1, nrow=2, labels="AUTO")
ggsave(filename='accuracy_figure.png', plot=gsa, device=NULL, width = 9.25, height=4.5, dpi=500)
ggsave(filename='accuracy_figure.pdf', plot=gsa, device=NULL, width = 9.25, height=4.5, dpi=600)
