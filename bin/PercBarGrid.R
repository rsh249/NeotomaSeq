#!/usr/bin/env Rscript
#args <- commandArgs(trailingOnly=TRUE)
setwd("~/Desktop/packrats")

#Load libraries
library(ggplot2)
library(dplyr)
library(data.table)
library(reshape2)
library(RColorBrewer)

#Pull files into a list based on classification level e.g. "Superkingdom"
names <- Sys.glob("*Superkingdom")
taxa <- lapply(names, read.delim, sep = "\t")

#Pull in packrat master data table to use later
full_table = read.delim("~/Desktop/packrats/packrat_aDNA.csv", sep = ",",  na.strings = "null")

#Create an empty data frame to fill with individual classification data
dfinal <- data.frame(Group.1 = character(),
                     numReads=integer(), 
                     numUniqueReads=integer(),
                     Extraction = character(),
                     Age = integer(),
                     Total = integer(),
                     stringsAsFactors=FALSE) 

#For each file, match the extraction information to the name of the file
for(i in 1:length(taxa)){
  df = taxa[[i]]
  #Split input name on '/'
  folders <- strsplit(names[i], '/')
  #Split again on '.'
  file <- strsplit(as.character(folders[length(folders)]), '.', fixed = TRUE)
  #Pull the extraction name from the file name to be matched up in the master data table
  extraction_name = as.character(file[[1]][1])
  match_name = substr(extraction_name, 1, nchar(extraction_name)-1) 
  
  #Add up the number of unique reads by classification
  df <- aggregate(cbind(count = numUniqueReads) ~ Group.1, 
                  data = df, 
                  FUN = sum)
  #Pull in file of filter read counts
  total.count <- read.delim("~/Desktop/packrats/read_counts.tab", sep = "\t", header = FALSE)
  
  #Match the extraction name with the file name
  df$Extraction <- NA
  df <- df %>%
    mutate(Extraction = match_name)
  
  #Create a column for the total reads for that sample
  df$Total <- total.count$V2[match(df$Extraction, total.count$V1)]
  
  #Remove processed samples
  full_table <- full_table[!grepl("*Pr", full_table$Name),]
  selected_rows <- full_table[grep(match_name, full_table$Name), ]
  
  #Pull the age from the master data table
  df$Age <- NA
  df <- df %>%
    mutate(Age = selected_rows$Midden.age[1])
  
  #Remove "NA" rows
  df <- df[!is.na(df$Age), ]
  
  #add individual sample data to larger data table
  dfinal <- rbind(dfinal, df)
}

#Average extractions from the same year
dfinal <- dfinal %>%
  group_by(Extraction, Group.1) %>%
  mutate_if(is.numeric, mean)

#Put on percentage scale
for(i in 1:nrow(dfinal)) {
  dfinal[i,] <- dfinal[i,] %>%
    mutate(count = count/Total*100)
  
}

#Only include one set data for each sample
dfinal <- unique( dfinal[ , 1:ncol(dfinal) ] )

AVV <- subset(dfinal, Group.1 %in% c("Archaea", "Viroids", "Viruses"))
BE <- subset(dfinal, Group.1 %in% c("Bacteria", "Eukaryota"))

#Graph the resulting table as a stacked bar graph colored by kingdom
png(filename =  "PercBarPlot.png", height = 7, width = 11, units = "in", res = 400)
#Order the table by midden age
dfinal$Age = factor(c(as.character(dfinal$Age)),levels=c("345","2835", "3105", "3260","3545", "28460", "31760"))
ggplot(data = dfinal, aes(x = Age, y = count, fill = Group.1)) + 
  labs(x = "Age (Years)", y = "Percentage of Filtered Reads") + 
  scale_x_discrete(labels=c("345 (COR)","2835 (COR)", "3105 (COR)", "3260 (COR)","3545 (GC)", "28460 (COR)", "31760 (COR)")) +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c("Archaea" = "#FDB462", "Bacteria" = "#80B1D3", "Eukaryota" = "#FB8072", "Viroids" = "#CCEBC5", "Viruses" = "#BEBADA"), name = "Kingdom") +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14), legend.text=element_text(size=12))
dev.off()

#Graph the resulting table as a stacked bar graph colored by kingdom -- Archaea, Viruses, Viroids
png(filename =  "AVVBarPlot.png", height = 7, width = 11, units = "in", res = 400)
#Order the table by midden age
AVV$Age = factor(c(as.character(AVV$Age)),levels=c("345","2835", "3105", "3260","3545", "28460", "31760"))
AVVplot = ggplot(data = AVV, aes(x = Age, y = count, fill = Group.1)) + 
  labs(x = "Age (Years)", y = "Percentage of Filtered Reads") + 
  scale_x_discrete(labels=c("345 (COR)","2835 (COR)", "3105 (COR)", "3260 (COR)","3545 (GC)", "28460 (COR)", "31760 (COR)")) +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c("Archaea" = "#FDB462", "Viroids" = "#CCEBC5", "Viruses" = "#BEBADA"), name = "Kingdom") +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=10), legend.text=element_text(size=8))
AVVplot
dev.off()

#Graph the resulting table as a stacked bar graph colored by kingdom -- Bacteria and Eukaryota
png(filename =  "BEBarPlot.png", height = 7, width = 11, units = "in", res = 400)
#Order the table by midden age
BE$Age = factor(c(as.character(BE$Age)),levels=c("345","2835", "3105", "3260","3545", "28460", "31760"))
BEplot = ggplot(data = BE, aes(x = Age, y = count, fill = Group.1)) + 
  labs(x = "Age (Years)", y = "Percentage of Filtered Reads") +
  scale_x_discrete(labels=c("345 (COR)","2835 (COR)", "3105 (COR)", "3260 (COR)","3545 (GC)", "28460 (COR)", "31760 (COR)")) +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c("Bacteria" = "#80B1D3", "Eukaryota" = "#FB8072"), name = "Kingdom") +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=10), legend.text=element_text(size=8))
BEplot
dev.off()
