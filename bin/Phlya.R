#Load libraries
library(ggplot2)
library(dplyr)
library(data.table)
library(reshape2)
library(RColorBrewer)

#Grab the TaxaResults file for each of the samples and look at the eukaryotic reads only
taxa <- read.delim("~/Desktop/packrats/PPCUn3kA_TAX/Taxa.ResultsAfter", sep = "\t")
eutaxa <- taxa[grepl("Eukaryota", taxa$superkingdom),]

#Count the number of reads for each phylum
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
  mutate(Extraction = "PPC529-Un")

#Assign the midden age
eutaxa$Age <- NA
eutaxa <- eutaxa %>%
  mutate(Age = as.numeric(3260))

#Repeat for all of the files and bind all them together
taxa2 <- read.delim("~/Desktop/packrats/PPCUn3kB_TAX/Taxa.ResultsAfter", sep = "\t")
eutaxa2 <- taxa2[grepl("Eukaryota", taxa2$superkingdom),]

eutaxa2 <- aggregate(cbind(count = numUniqueReads) ~ phylum, 
                    data = eutaxa2, 
                    FUN = sum)

total.count <- sum(eutaxa2$count)

eutaxa2$Total <- NA
eutaxa2 <- eutaxa2 %>%
  mutate(Total = total.count)

eutaxa2$Extraction <- NA
eutaxa2 <- eutaxa2 %>%
  mutate(Extraction = "PPC529-Un")

eutaxa2$Age <- NA
eutaxa2 <- eutaxa2 %>%
  mutate(Age = as.numeric(3260))

eutaxa <- rbind(eutaxa, eutaxa2)


taxa3 <- read.delim("~/Desktop/packrats/TSUn28kA_TAX/Taxa.ResultsAfter", sep = "\t")
eutaxa3 <- taxa3[grepl("Eukaryota", taxa3$superkingdom),]

eutaxa3 <- aggregate(cbind(count = as.numeric(numUniqueReads)) ~ phylum, 
                     data = eutaxa3, 
                     FUN = sum)

total.count <- sum(eutaxa3$count)

eutaxa3$Total <- NA
eutaxa3 <- eutaxa3 %>%
  mutate(Total = total.count)

eutaxa3$Extraction <- NA
eutaxa3 <- eutaxa3 %>%
  mutate(Extraction = "TS211F-Un")

eutaxa3$Age <- NA
eutaxa3 <- eutaxa3 %>%
  mutate(Age = as.numeric(28460))

eutaxa <- rbind(eutaxa, eutaxa3)


taxa4 <- read.delim("~/Desktop/packrats/FRT504", sep = "\t")
eutaxa4 <- taxa4[grepl("Eukaryota", taxa4$superkingdom),]

eutaxa4 <- aggregate(cbind(count = as.numeric(numUniqueReads)) ~ phylum, 
                     data = eutaxa4, 
                     FUN = sum)

total.count <- sum(eutaxa4$count)

eutaxa4$Total <- NA
eutaxa4 <- eutaxa4 %>%
  mutate(Total = total.count)

eutaxa4$Extraction <- NA
eutaxa4 <- eutaxa4 %>%
  mutate(Extraction = "FRT504-Un")

eutaxa4$Age <- NA
eutaxa4 <- eutaxa4 %>%
  mutate(Age = as.numeric(345))

eutaxa <- rbind(eutaxa, eutaxa4)


taxa5 <- read.delim("~/Desktop/packrats/FRT511A", sep = "\t")
eutaxa5 <- taxa5[grepl("Eukaryota", taxa5$superkingdom),]

eutaxa5 <- aggregate(cbind(count = as.numeric(numUniqueReads)) ~ phylum, 
                     data = eutaxa5, 
                     FUN = sum)

total.count <- sum(eutaxa5$count)

eutaxa5$Total <- NA
eutaxa5 <- eutaxa5 %>%
  mutate(Total = total.count)

eutaxa5$Extraction <- NA
eutaxa5 <- eutaxa5 %>%
  mutate(Extraction = "FRT511A-Un")

eutaxa5$Age <- NA
eutaxa5 <- eutaxa5 %>%
  mutate(Age = as.numeric(2835))

eutaxa <- rbind(eutaxa, eutaxa5)


taxa6 <- read.delim("~/Desktop/packrats/PPC524", sep = "\t")
eutaxa6 <- taxa6[grepl("Eukaryota", taxa6$superkingdom),]

eutaxa6 <- aggregate(cbind(count = as.numeric(numUniqueReads)) ~ phylum, 
                     data = eutaxa6, 
                     FUN = sum)

total.count <- sum(eutaxa6$count)

eutaxa6$Total <- NA
eutaxa6 <- eutaxa6 %>%
  mutate(Total = total.count)

eutaxa6$Extraction <- NA
eutaxa6 <- eutaxa6 %>%
  mutate(Extraction = "PPC524-Un")

eutaxa6$Age <- NA
eutaxa6 <- eutaxa6 %>%
  mutate(Age = as.numeric(3105))

eutaxa <- rbind(eutaxa, eutaxa6)


taxa7 <- read.delim("~/Desktop/packrats/TS211A", sep = "\t")
eutaxa7 <- taxa7[grepl("Eukaryota", taxa7$superkingdom),]

eutaxa7 <- aggregate(cbind(count = as.numeric(numUniqueReads)) ~ phylum, 
                     data = eutaxa7, 
                     FUN = sum)

total.count <- sum(eutaxa7$count)

eutaxa7$Total <- NA
eutaxa7 <- eutaxa7 %>%
  mutate(Total = total.count)

eutaxa7$Extraction <- NA
eutaxa7 <- eutaxa7 %>%
  mutate(Extraction = "TS211A-Un")

eutaxa7$Age <- NA
eutaxa7 <- eutaxa7 %>%
  mutate(Age = as.numeric(31760))

eutaxa <- rbind(eutaxa, eutaxa7)


#Average extractions from the same year
avg <- eutaxa %>%
  group_by(Extraction,phylum) %>%
  mutate_if(is.numeric, mean)

for(i in 1:nrow(avg)) {
  avg[i,] <- avg[i,] %>%
    mutate(count = count/Total*100)
  
}

avg <- unique( avg[ , 1:ncol(avg) ] )
avg = avg[order(-avg$count),]
avg = avg[order(avg$Extraction),]

newavg <- do.call("rbind", list(avg[1:5,], avg[44:48,], avg[81:85,], avg[126:130,], avg[176:180,], avg[223:226,]))

newavg$Age = factor(c(newavg$Age),levels=c("345","2835","3105","3260", "28460", "31760"))

#Graph the resulting table as a stacked bar graph colored by kingdom on a logarithmic scale
png(filename =  "PhylaBarPlot.png", height = 7, width = 11, units = "in", res = 400)

#dfinal$Group.1 <- factor(dfinal$Group.1, levels = unique(dfinal$Group.1[order(-as.numeric(as.character(dfinal$Age)))]))
phylaplot = ggplot(data = newavg, aes(x = Age, y = count, fill = phylum)) + 
  labs(x = "Age (Years)", y = "Percentage of Unique Eukaryotic Reads") + 
  scale_x_discrete(labels=c("345 (COR)","2835 (COR)", "3105 (COR)", "3260 (COR)", "28460 (COR)", "31760 (COR)")) +
  geom_bar(stat = "identity") + 
  #scale_fill_discrete(name = "Phylum") +
  scale_fill_brewer(direction = -1, palette="Set3", name = "Phylum") + 
  theme(axis.text=element_text(size=8), axis.title=element_text(size=10), legend.text=element_text(size=8))
phylaplot
dev.off()

