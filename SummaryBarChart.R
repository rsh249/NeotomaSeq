#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

library(ggplot2)
library(dplyr)
library(data.table)
library(reshape2)

names <- Sys.glob(args[[1]])
print(names)
taxa <- lapply(names, read.delim, sep = "\t")

full_table = read.delim(args[2], sep = ",",  na.strings = "null")
dfinal <- data.frame(Group.1 = character(),
                 numReads=integer(), 
                 numUniqueReads=integer(),
                 Extraction = character(),
                 Age = integer(),
                 stringsAsFactors=FALSE) 
for(i in 1:length(taxa)){
  df = taxa[[i]]
  
  folders <- strsplit(names[i], '/')
  file <- strsplit(as.character(folders[length(folders)]), '.', fixed = TRUE)
  
  extraction_name = as.character(file[[1]][1])
  match_name = substr(extraction_name, 1, nchar(extraction_name)-1) 
  df$Extraction <- NA
  df <- df %>%
    mutate(Extraction = match_name)
  #Remove processed samples
  full_table <- full_table[!grepl("*Pr", full_table$Name),]
  selected_rows <- full_table[grep(match_name, full_table$Name), ]
  df$Age <- NA
  df <- df %>%
    mutate(Age = selected_rows$Midden.age[1])
  #Remove "NA" rows
  df <- df[!is.na(df$Age), ]
  dfinal <- rbind(dfinal, df)
}

#Average extractions from the same year
dfinal <- dfinal %>%
  group_by(Extraction,Group.1) %>%
  mutate_if(is.numeric, mean)

print(dfinal)

#Graph the resulting table as a stacked bar graph colored by kingdom on a logarithmic scale
png(filename =  "SummaryBarPlot.png", height = 7, width = 11, units = "in", res = 600)
#Order the table
dfinal = dfinal[order(-dfinal$Age),]
#dfinal$Group.1 <- factor(dfinal$Group.1, levels = unique(dfinal$Group.1[order(-as.numeric(as.character(dfinal$Age)))]))
ggplot(data = dfinal, aes(x = as.numeric(as.character(Age)), y = numUniqueReads, fill = Group.1)) + labs(x = "Age (Years)", y = "Number of Unique Reads (log scale)") + geom_bar(stat = "identity") + scale_y_continuous(trans='log10')
dev.off()
