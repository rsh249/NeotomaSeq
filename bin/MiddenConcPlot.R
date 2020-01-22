library(ggplot2)
library(dplyr)
library(stringr)
library(reshape2)

#Read in table of packrat information
p.data = read.delim("Files2Run/packrat_aDNA_notes.csv", sep = ",")
#Make the bioanalyzer column numeric
p.data$Bioanalyzer.conc.ng.ul <- as.numeric(as.character(p.data$Bioanalyzer.conc.ng.ul))
#Remove crazy outlier
p.data <- p.data[-19,]
op.data = p.data
#Order it by midder Age
op.data$Midden.age = as.numeric(op.data$Midden.age)
op.data = p.data[order(p.data$Midden.age),]
op.data$samName = as.character(op.data$Name)
op.data$Site = stringr::word(op.data$Midden.loc, 1, 2) %>% 
  stringr::str_replace_all(" ", "_") %>% 
  stringr::str_replace_all("Flaming_Rock", "City of Rocks") %>% 
  stringr::str_replace_all("Pinnacle_Pass", "City of Rocks") %>%
  stringr::str_replace_all("Twin_Sisters", "City of Rocks") %>%  
  tidyr::replace_na('Blank')

#Plot results
(concfig = ggplot(data = op.data, aes(x = reorder(paste(X, samName, sep = '_'), Midden.age), y = as.numeric(as.character(Qubit.conc.ng.ul)), fill = Site)) + 
  geom_col( width=0.7) + 
  coord_flip() +
  geom_text(data = op.data, aes(label=Midden.age), angle=0, vjust=0.5, hjust=0, size=2) +
  theme(axis.text.x = element_text(angle=45,hjust=1, size=6), 
        axis.text.y = element_text(angle = 45, size =6, hjust=1),
        legend.text=element_text(size=6), 
        legend.title=element_text(size=8),
        axis.title=element_text(size=8),
        panel.background = element_rect(fill = 'white', colour = 'white'),  
        panel.grid.major = element_line(size = 0.25, linetype = 'solid',colour = "gray"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "gray"),
        legend.position = 'bottom') + 
  xlab("Extraction Label (by Midden Age)") + 
  ylab(bquote("aDNA Concentration" ~ (ng/mu*L)))
)

ggsave(filename = "ConcPlot.pdf", plot=concfig, device=NULL, height = 9.25, width = 4.5, dpi = 500)


