library(ggplot2)
library(dplyr)
library(stringr)
library(reshape2)

#Read in table of packrat information
p.data = read.delim("Desktop/packrats/packrat_aDNA.csv", sep = ",")
#Make the bioanalyzer column numeric
p.data$Bioanalyzer.conc.ng.ul <- as.numeric(as.character(p.data$Bioanalyzer.conc.ng.ul))
#Remove crazy outlier
p.data <- p.data[-19,]
#Order it by midder Age
op.data = p.data[order(p.data$Midden.age),]
#Plot results
png(filename = "ConcPlot.png", height = 7, width = 11, units = "in", res = 400)
ggplot(data = op.data, aes(x = reorder(Name, Midden.age), y = Bioanalyzer.conc.ng.ul, fill = Site)) + geom_bar(stat = "identity") + 
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5, size=8), legend.text=element_text(size=12), 
        panel.background = element_rect(fill = 'white', colour = 'white'),  
        panel.grid.major = element_line(size = 0.25, linetype = 'solid',colour = "gray"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "gray")) + 
  xlab("Extraction Label (by Midden Age)") + 
  ylab(bquote("aDNA Concentration" ~ (ng/mu*L)))
dev.off()
