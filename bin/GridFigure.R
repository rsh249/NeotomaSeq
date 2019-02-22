library(ggplot2)
library(gridExtra)

source("./PercBarGrid.R")
source("./Phlya.R")
source("./Families.R")

png(filename =  "GridPlot.png", height = 9, width = 15, units = "in", res = 400)
grid.arrange(BEplot, AVVplot, phylaplot, famplot, nrow = 2)
dev.off()
