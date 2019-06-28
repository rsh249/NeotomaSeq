# read data for:

# Macrofossils
macro = read.table("./COR_macrofossils.csv", sep=',', header =T, stringsAsFactors = F)

colnames(macro) = macro[1,]
macro=macro[-1,]
head(macro)

#compile families by site
n=1
genera = list()
sites = vector()
for(i in 5:ncol(macro)){
  names = vector()
  for(z in macro$Name[which(macro[,i]>0)]){
    names = append(names, unique(strsplit(z, ' ')[[1]][1]))
  }
  genera[[n]] = names
  sites[n] = colnames(macro)[i]
  n=n+1
}
print(genera)

## Read the aDNA genus files
dnagenfiles = list.files('~/NeotomaSeq/results', recursive=TRUE, pattern='Taxa.ResultsAfter', full.names = T)
dnagenera = list()
dnasites = vector()
for(i in dnagenfiles){
  nom = strsplit(i, '/Taxa.ResultsAfter')[[1]][1]
  nom2 = strsplit(nom, '[/]')[[1]]
  name = nom2[length(nom2)]
  dnagen = data.table::fread(i, sep = '\t')
  #get plant hits ONLY
  dnasub = dnagen[dnagen$phylum == 'Streptophyta',]
  sumreads = sum(dnagen$numUniqueReads)
  aggr = aggregate(dnasub$numUniqueReads, by=list(dnasub$genus), FUN=sum )
  subs = aggr[aggr$x>(0.001*sumreads),]
  dnasites[name] = name
  dnagenera[[name]] = as.character(subs$Group.1)
}

## Read the aDNA family files
dnafamfiles = list.files('~/NeotomaSeq/results', recursive=TRUE, pattern='Taxa.ResultsAfter', full.names = T)
dnafamily = list()
dnafamsites = vector()
for(i in dnafamfiles){
  nom = strsplit(i, '/Taxa.ResultsAfter')[[1]][1]
  nom2 = strsplit(nom, '[/]')[[1]]
  name = nom2[length(nom2)]
  dnafam = data.table::fread(i, sep = '\t')
  #get plant hits ONLY
  dnasub = dnafam[dnafam$phylum == 'Streptophyta',]
  sumreads = sum(dnagen$numUniqueReads)
  aggr = aggregate(dnasub$numUniqueReads, by=list(dnasub$family), FUN=sum )
  subs = aggr[aggr$x>(0.001*sumreads),]
  dnafamsites[name] = name
  dnafamily[[name]] = as.character(subs$Group.1)
}

#Amplicon sequencing (ITS2/rbcL)
ampl = list.files('compiled_results/', full.names = T)
ampl.nums = vector()
ampl.barcode=vector()
ampl.gen = list()
for(a in ampl){
  path = strsplit(a, '[/]')[[1]]
  name = path[3]
  nsplit = strsplit(name, '[_]')[[1]]
  num=nsplit[3]
  gene=strsplit(nsplit[4], '[-.]')[[1]][2]
  hits = data.table::fread(a)
  #print(hits)
  print(gene)
  firstgen=vector()
  for(k in 1:nrow(hits)){
    firsthit = strsplit(hits$names[k], '[;]')[[1]][1]
    firstgen[k] = strsplit(firsthit, '[ ]')[[1]][1]
  }
  ampl.nums[name] = num
  ampl.gen[[name]] = c(firstgen)
  ampl.barcode[name] = gene
}

#Merge by sample into a single table
nummap = list('PPC502-Un', 'TS211A-Un', 'FRT531-Un', 'FRT511A-Un', 'TS564-Un', 'GC100B-Un', 'Blank')
names(nummap) = c('01', '05', '09', '13', '17', '21', '0N')
ampl.sites = unlist(nummap[ampl.nums])
# PPC502-Un = 01 **OLDEST SAMPLE @ 47500
# TS211A-Un = 05
# FRT531-Un = 09
# FRT511A-Un = 13
# TS564-Un = 17
# GC100B-Un = 21
# Blank = 0N

allfoss=unique(unlist(genera))
nr = length(allfoss)
coll = matrix(nrow=nr, ncol = 18)
xx = 1
columns = vector()
for(nn in 1:length(dnasites)){
  q=dnasites[[nn]]
  
  presite = strsplit(q, '[-]')[[1]][1]
  r = grep(presite, sites)
  if(length(r)==0){next}
  columns[[xx]] = paste(sites[[r]], '.foss', sep='')
  coll[which(allfoss %in% genera[[r]]),xx] = 1
  
  
  columns[[xx]] = paste(dnasites[[nn]], '.dnagen', sep ='')
  coll[which(allfoss %in% dnafamily[[nn]]),] = 1
  xx=xx+1
  
  
  s=grep(presite,ampl.sites)
  #merge gene and site name to amplicon identifier
  paste(ampl.barcode[s], sep = '_')
  columns[[xx]] = amplicon.id
  #add amplicon hits to table
  coll[which(allfoss %in% ampl.gen[s]),]=1
  #add amplicon family hits to table
  #taxonomizr
  #sub.fams = translate.genus.to.family
  coll[which(allfoss %in% sub.fams),]=1
  
  
}


#plot using pheatmap or ggplot2 tile




