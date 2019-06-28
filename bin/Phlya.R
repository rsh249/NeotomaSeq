#Load libraries
library(ggsci)
library(ggplot2)
library(dplyr)
library(data.table)
library(reshape2)
library(RColorBrewer)
root = "~/NeotomaSeq"
files = list.files(root, recursive=T, full.names = T, pattern="*.ResultsAfter")

full_table = read.delim("./Files2Run/packrat_aDNA.csv", sep = ",",  na.strings = "null")



eutaxadf = data.frame(family = character(),
                      count=numeric(), 
                      Total=numeric(),
                      Extraction = character(),
                      Age = numeric(),
                      stringsAsFactors=FALSE)
for(i in files){
  #parse name
  lastname=name
  n1 = strsplit(i, '[/]')[[1]][6]
  n2 = strsplit(n1, '[-]')
  name = n2[[1]][1]
  
  if(name == "GC100B"){next} #Make another figure to compare MX and US sites
  
  #grep name against full_table to get age data
  age = unique(full_table[grep(name, full_table$Name), 'Midden.age'])
  name = paste(name, n2[[1]][2], sep = '-') #put name identifier back together
  print(name)
  
  #Grab the TaxaResults file for each of the samples and look at the eukaryotic streptophyla reads only
  taxa <- data.table::fread(i, sep = "\t")
  eutaxa <- taxa[grepl("Eukaryota", taxa$superkingdom),]
  #eutaxa <- taxa[grepl("Streptophyta", taxa$phylum),]
  
  #Count the number of reads for each family
  eutaxa <- aggregate(cbind(count = numUniqueReads) ~ phylum, 
                      data = eutaxa, 
                      FUN = sum)
  #Sum the total eukarotic reads
  total.count <- sum(eutaxa$count)
  
  #Put the total count of eukaryotic reads in the eutaxa table
  eutaxa$Total <- NA
  eutaxa <- eutaxa %>%
    mutate(Total = total.count)
  
  #Label each extraction based on the file name
  eutaxa$Extraction <- NA
  eutaxa <- eutaxa %>%
    mutate(Extraction = n1)
  #Assign the midden age
  eutaxa$Age <- NA
  eutaxa <- eutaxa %>%
    mutate(Age = as.numeric(as.character(age)))
  eutaxadf=rbind(eutaxadf, eutaxa)
}


#Trim to top 15 families across samples

trimtaxa = unique(eutaxadf[order(-(eutaxadf$count/eutaxadf$Total)),'phylum'][1:100])[1:8]
ttaxadf = eutaxadf[which(eutaxadf$phylum %in% trimtaxa),]

ttaxadf$plotage = cut(ttaxadf$Age, breaks=c(0, 1000, 5000, 10000, 15000, 20000, 25000, Inf), 
                      labels=c('<1ka', '1-5ka', '5-10ka', '10-15ka', '15-20ka', '20-25ka', '>25ka'))

phylaplot = ggplot(data = ttaxadf) + 
  geom_col(aes(x = reorder(Extraction, Age), y = count/Total*100, fill = phylum)) + 
  scale_color_npg()+ scale_fill_npg() + 
  theme_linedraw() + 
  theme(
    axis.text = element_text(size = 8),
    axis.text.x  = element_text(angle = 90, size =0),
    axis.title = element_text(size = 8),
    legend.text = element_text(size = 8)
  ) + labs(x = "", y = "% Unique Eukaryote Reads", fill='Phylum') +
  facet_grid(.~plotage, scales = 'free_x', space='free') 
phylaplot
#dev.off()

