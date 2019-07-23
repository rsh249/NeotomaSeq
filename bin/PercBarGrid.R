#!/usr/bin/env Rscript
#args <- commandArgs(trailingOnly=TRUE)
#setwd("~/Desktop/packrats")

#Load libraries
library(ggplot2)
library(ggsci)
library(dplyr)
library(data.table)
library(reshape2)
library(RColorBrewer)
root='~/NeotomaSeq'

#Pull files into a list based on classification level e.g. "Superkingdom"
names <- Sys.glob(paste(root, "/results/*/*Superkingdom", sep = ''))
taxa <- lapply(names, read.delim, sep = "\t")

#Pull in packrat master data table to use later
full_table = read.delim("./Files2Run/packrat_aDNA.csv", sep = ",",  na.strings = "null")

#Pull in file of filter read counts
total.count <- read.delim("./Files2Run/read_counts.tab", sep = "\t", header = TRUE)

#Create an empty data frame to fill with individual classification data
dfinal <- data.frame(Group.1 = character(),
                     numReads=numeric(), 
                     numUniqueReads=numeric(),
                     Extraction = character(),
                     Age = numeric(),
                     Total = numeric(),
                     stringsAsFactors=FALSE) 

#For each file, match the extraction information to the name of the file
for(i in 1:length(taxa)){
  data = taxa[[i]]
  n1 = strsplit(names[i], '[/]')[[1]][6]
  n2 = strsplit(n1, '[-]')
  match_name = n2[[1]][1]
  full_name = paste(match_name, n2[[1]][2], sep='-')
  print(n1)
  print(full_name)


  if(match_name == "GC100B"){next} #Make another figure to compare MX and US sites
  
  #grep name against full_table to get age data
  age = unique(full_table[grep(match_name, full_table$Name), 'Midden.age'])
  
  #Add up the number of unique reads by classification
  df <- aggregate(cbind(count = numUniqueReads) ~ Group.1, 
                  data = data, 
                  FUN = sum)

  #Match the extraction name with the file name
  df[,'Extraction'] = rep(n1, nrow(df))
  
  #Create a column for the total reads for that sample
  #df$Total[i] = total.count[grep(match_name, as.character(total.count$V1)),1]
  df[,'Total'] = rep(sum(data$numUniqueReads, na.rm=T), nrow(df))
  
  #Remove processed samples
 # full_table <- full_table[!grepl("*Pr", full_table$Name),]
 
 selected_rows <- full_table[grep(match_name, full_table$Name), ]

  #Pull the age from the master data table
  df$Age <- as.numeric(as.character(selected_rows$Midden.age[1]))
  
  #Remove "NA" rows
  df <- df[!is.na(df$Age), ]
  
  #add individual sample data to larger data table
  dfinal <- rbind(dfinal, df)
 # print(dfinal)
}


#Average extractions from the same year
#dfinal <- dfinal %>%
#  group_by(Extraction, Group.1) %>%
#  mutate_if(is.numeric, mean)

#Put on percentage scale
for(i in 1:nrow(dfinal)) {
  dfinal[i,'count'] <- dfinal$count[i]/as.numeric(dfinal$Total[i])*100

}

#Only include one set data for each sample
dfinal <- unique( dfinal[ , 1:ncol(dfinal) ] )

#Graph the resulting table as a stacked bar graph colored by kingdom
#png(filename =  "PercBarPlot.png", height = 7, width = 11, units = "in", res = 400)
#Order the table by midden age
dfinal$plotage = cut(dfinal$Age, breaks=c(0, 1000, 5000, 10000, 15000, 20000, 25000, Inf), 
                      labels=c('<1ka', '1-5ka', '5-10ka', '10-15ka', '15-20ka', '20-25ka', '>25ka'))
#dfinal$Age = factor(c(as.character(dfinal$Age)),levels=c("345","2835", "3105", "3260","3545", "28460", "31760"))

AVV <- subset(dfinal, Group.1 %in% c("Archaea", "Viroids", "Viruses"))
BE <- subset(dfinal, Group.1 %in% c("Bacteria", "Eukaryota"))

#Just plot bacteria and Eukaryotes
BEplot = ggplot(data = BE) +
  geom_col( aes(
    x = reorder(Extraction,Age),
    y = count,
    fill = Group.1
  )) + scale_color_npg()+ scale_fill_npg() + 
  theme_linedraw() + 
  theme(
    axis.text = element_text(size = 8),
    axis.text.x  = element_text(angle = 90, size =0),
    axis.title = element_text(size = 8),
    legend.text = element_text(size = 8)
  ) + 
  labs(x = "", y = "% Classified Reads", fill='Kingdom') +
  facet_grid(.~plotage, scales = 'free_x', space='free') + ylim(c(0,100))

BEplot


#Graph the resulting table as a stacked bar graph colored by kingdom -- Archaea, Viruses, Viroids
#png(filename =  "AVVBarPlot.png", height = 7, width = 11, units = "in", res = 400)
#Order the table by midden age
AVV=AVV[which(AVV$Group.1!='Viroids'),]
AVVplot = ggplot(data = AVV) +
  geom_col( aes(
    x = reorder(Extraction, Age),
    y = count,
    fill = Group.1
  )) + scale_color_npg()+ scale_fill_npg() + 
  theme_linedraw() + 
  theme(
    axis.text = element_text(size = 8),
    axis.text.x  = element_text(angle = 90, size =0),
    axis.title = element_text(size = 8),
    legend.text = element_text(size = 8)
  ) + 
  labs(x = "", y = "% Classified Reads", fill='Kingdom') +
  facet_grid(.~plotage, scales = 'free_x', space='free') 

AVVplot
#dev.off()


