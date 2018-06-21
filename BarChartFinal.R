#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

library(ggplot2)
library(dplyr)

#Pull in a Tax.Results file and sum the number of reads for each kingdom
df = read.delim(args[[1]], sep = "\t")
df <- aggregate(cbind(count = Num_Unique_Reads) ~ Kingdoms, 
                data = df, 
                FUN = sum)
#Label the extraction
df$Extraction<-NA
df <- df %>%
  mutate(Extraction = basename(args[[1]]))

#Pull in a Tax.Results file and sum the number of reads for each kingdom
df2  = read.delim(args[[2]], sep = "\t")
df2 <- aggregate(cbind(count = Num_Unique_Reads) ~ Kingdoms, 
                 data = df2, 
                 FUN = sum)
#Label the extraction
df2$Extraction<-NA
df2 <- df2 %>%
  mutate(Extraction = basename(args[[2]]))

#Bind the tables together
df <- rbind(df, df2)

#Pull in a Tax.Results file and sum the number of reads for each kingdom
df3  = read.delim(args[[3]], sep = "\t")
df3 <- aggregate(cbind(count = Num_Unique_Reads) ~ Kingdoms, 
                 data = df3, 
                 FUN = sum)
#Label the extraction
df3$Extraction<-NA
df3 <- df3 %>%
  mutate(Extraction = basename(args[[3]]))

#Bind the tables together
df <- rbind(df, df3)

#Pull in a Tax.Results file and sum the number of reads for each kingdom
df4  = read.delim(args[[4]], sep = "\t")
df4 <- aggregate(cbind(count = Num_Unique_Reads) ~ Kingdoms, 
                 data = df4, 
                 FUN = sum)
#Label the extraction
df4$Extraction<-NA
df4 <- df4 %>%
  mutate(Extraction = basename(args[[4]]))

#Bind the tables together
df <- rbind(df, df4)

#Order the table
df <- df[order(-df$count),]

#Graph the resulting table as a stacked bar graph colored by kingdom on a logarithmic scale
ggplot(data = df, aes(x = Extraction, y = count, fill = Kingdoms)) + geom_bar(stat = "identity") + scale_y_continuous(trans='log10')

