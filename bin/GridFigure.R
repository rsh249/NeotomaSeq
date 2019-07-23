library(ggplot2)
library(gridExtra)

source("./PercBarGrid.R")
source("./Phlya.R")
source("./Families.R")
source("./gen_compare_fossil_DNA.R")
source('./accuracy_map.R')

png(filename =  "GridPlot.png", height = 11, width = 8, units = "in", res = 400)
grid.arrange(BEplot, AVVplot, phylaplot, famplot, nrow = 4, heights=c(1,1,1,1.5))
dev.off()



pp1 <- ggplotGrob(BEplot)

pp2 <- ggplotGrob(AVVplot)

pp3 <- ggplotGrob(phylaplot)

pp4 <- ggplotGrob(famplot)



library(cowplot)
png(filename = 'GridPlot3.png', height = 9.5, width = 5.75, units = 'in', res=400)
plot_grid(pp1, pp2, pp3, pp4, ncol=1, align = 'v', rel_heights=c(0.85,0.85,1.3,1.45), labels="AUTO")
dev.off()
