# read data for:

# Macrofossils
macro = read.table("./COR_macrofossils.csv", sep=',', header =T, stringsAsFactors = F)

colnames(macro) = macro[1,]
macro=macro[-1,]
head(macro)

#compile families by site
n=1
fams = list()
sites = vector()
for(i in 5:ncol(macro)){
  fams[[n]] = (unique(macro$Family[which(macro[,i]>0)]))
  sites[n] = colnames(macro)[i]
  n=n+1
}

## Read the aDNA family files
dirs = list.dirs('~/NeotomaSeq/results')
# files = list.files(dirs, pattern='.Family', full.names = T)
# eutaxadf = data.frame(genus = character(),
#                       count=numeric(), 
#                       Total=numeric(),
#                       Extraction = character(),
#                       Age = numeric(),
#                       stringsAsFactors=FALSE)
# name=''
# for(i in files){
#   #parse name
#   lastname=name
#   n1 = strsplit(i, '[/]')[[1]][6]
#   n2 = strsplit(n1, '[-]')
#   name = n2[[1]][1]
#   
#   if(name == "GC100B"){next} #Make another figure to compare MX and US sites
#   
#   #grep name against full_table to get age data
#   age = unique(full_table[grep(name, full_table$Name), 'Midden.age'])
#   name = paste(name, n2[[1]][2], sep = '-') #put name identifier back together
#   print(name)
#   
#   #Grab the TaxaResults file for each of the samples and look at the eukaryotic streptophyla reads only
#   taxa <- data.table::fread(i, sep = "\t")
#   eutaxa = taxa
#   #eutaxa <- taxa[grepl("Eukaryota", taxa$superkingdom),]
#   #eutaxa <- taxa[grepl("Streptophyta", taxa$phylum),]
#   
#   #Count the number of reads for each family
#   #eutaxa <- aggregate(cbind(count = numUniqueReads) ~ Group.1, 
#                      # data = eutaxa, 
#                      # FUN = sum)
#   #Sum the total eukarotic reads
#   total.count <- sum(eutaxa$numUniqueReads)
#   
#   #Put the total count of eukaryotic reads in the eutaxa table
#   eutaxa$Total <- NA
#   eutaxa <- eutaxa %>%
#     mutate(Total = total.count)
#   
#   #Label each extraction based on the file name
#   eutaxa$Extraction <- NA
#   eutaxa <- eutaxa %>%
#     mutate(Extraction = n1)
#   
#   #Assign the midden age
#   eutaxa$Age <- NA
#   eutaxa <- eutaxa %>%
#     mutate(Age = as.numeric(as.character(age)))
#   eutaxadf=rbind(eutaxadf, eutaxa)
# }
# 
# eutaxadf = eutaxadf[order(-eutaxadf[,'numUniqueReads']),]



## Read the aDNA family files
dnafamfiles = list.files('~/NeotomaSeq/results', recursive=TRUE, pattern='Taxa.ResultsAfter', full.names = T)
dnafams = list()
dnafamsites = vector()
for(i in dnafamfiles){
  nom = strsplit(i, '/Taxa.ResultsAfter')[[1]][1]
  nom2 = strsplit(nom, '[/]')[[1]]
  name = nom2[length(nom2)]
  dnafam = data.table::fread(i, sep = '\t')
  #get plant hits ONLY
  dnasub = dnafam[dnafam$phylum == 'Streptophyta',]
  sumreads = sum(dnasub$numUniqueReads)
  aggr = aggregate(dnasub$numUniqueReads, by=list(dnasub$family), FUN=sum )
  subs = aggr[aggr$x>(0.01*sumreads),]
  dnafamsites[name] = name
  dnafams[[name]] = as.character(subs$Group.1)
}

#calc accuracy
dna2macro = matrix(ncol=10, nrow=5)
columns=vector()
zz=1

for(i in 1:length(sites)) {
  b = grep(sites[i], dnafamsites)
  if(length(b)==1) {
    #samplename = strsplit(dnasites[[b]], '/')[[1]][6]
    samplename=dnafamsites[[b]]
    columns[zz] = samplename
    summ = sum(fams[[i]] %in% dnafams[[b]])
    dna2macro[1,zz] = summ
    dna2macro[2,zz] = summ/length(fams[[i]])
    dna2macro[3,zz] = summ/length(dnafams[[b]])
    dna2macro[4,zz] = length(fams[[i]])
    dna2macro[5,zz] = length(dnafams[[b]])
    zz=zz+1
  } else {
    for(n in b){
      samplename = dnafamsites[[n]]
      columns[zz] = samplename
      summ = sum(fams[[i]] %in% dnafams[[n]])
      dna2macro[1,zz] = summ
      dna2macro[2,zz] = summ/length(fams[[i]])
      dna2macro[3,zz] = summ/length(dnafams[[n]])
      dna2macro[4,zz] = length(fams[[i]])
      dna2macro[5,zz] = length(dnafams[[n]])
      zz=zz+1
    }
  }
  
}

dna2macro=as.data.frame(dna2macro)
colnames(dna2macro) = columns
rownames(dna2macro) = c('intersect', 
                        'TPR_macro', 
                        'TPR_DNA', 
                        'totalmacro', 
                        'totalDNA')
dna2macro=as.data.frame(t(dna2macro))
dna2macro = cbind(columns, dna2macro)

p1 = ggplot(dna2macro) + geom_col(aes(x = columns, y=TPR_DNA)) + ylim(c(0,1))
p2 = ggplot(dna2macro) + geom_col(aes(x = columns, y=TPR_macro)) + ylim(c(0,1))

library(gridExtra)
ggsave('TPR.png', 
       grid.arrange(p1, p2))







