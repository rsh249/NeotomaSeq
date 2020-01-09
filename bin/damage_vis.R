

#Load libraries
library(ggsci)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(data.table)
library(reshape2)
library(RColorBrewer)
library(vroom)

full_table = read.delim("./Files2Run/packrat_aDNA.csv", sep = ",",  na.strings = "null")



# Collect 5' C->T transition rates from
# $damage_results/$sample/4pCtoT_freq.txt

samples = list.files('~/NeotomaSeq/damage_results', full.names = T)


# for (i in 1:length(samples)) {
#   s = samples[i]
# 
#   
#   r = read.table(paste(s, "/results_Juniperus.sort/5pCtoT_freq.txt",  sep=''), header =T)
#   rf = read.table(paste(s, "/results_Juniperus.pmds3filter/5pCtoT_freq.txt",  sep=''), header =T)
#   print(s)
#   #plot(r[,1], r[,2] - rf[,2], ylim=c(-0.2,0.2))
#   plot(r, ylim=c(0,0.2))
#   points(rf, col='navy', pch=20)
#   Sys.sleep(3)
#   #print(r)
# }


# collCT = data.frame(name=character(),
#                     FivepCT.mean=numeric(),
#                     FivepCT.sd=numeric(),
#                     MidCT.mean=numeric(),
#                     MidCT.sd=numeric(),
#                     Age=numeric())

##Assemble deamination tables to work in ggplot with names and ages attached

coll.ctp = data.frame(name=character(),
                      age=numeric(), 
                      z=numeric(),
                      CT=numeric(),
                      CA=numeric(),
                      CG=numeric(),
                      CC=numeric(),
                      GA=numeric(),
                      GT=numeric(),
                      GC=numeric(),
                      GG=numeric()
)


coll.readlength = data.frame(name=character(),
                             age=numeric(),
                             readlength=numeric(),
                             count=numeric()
                             )
for (i in 1:length(samples)) {
  if(grep('GC100', samples)==i){next}
  ### FIX NAMES TO MATCH OTHER FIGURES ###
  s = samples[i]
  s.split = strsplit(s, '/')
  name = strsplit(s.split[[1]][6], '-')[[1]][1]
  age = unique(full_table[grep(name, full_table$Name), 'Midden.age'])
  nt = s.split[[1]][6]

 print(s)
  ctp = read.table(paste(s, '/Plastid.deam', sep=''), header=T)
  ctp = cbind(rep(nt, nrow(ctp)), rep(age, nrow(ctp)), ctp)
  names(ctp) = c('name', 'age', 'z', 'CT', 'CA', 'CG', 'CC', 'GA', 'GT', 'GC', 'GG')
  

  readl = vroom(paste(s, '/read_length.txt', sep=''), col_types = 'n', num_threads=2, col_names=F)
  readl = readl %>% 
    mutate(a = cut(X1, breaks =seq(0,250,10), labels=seq(10,250,10))) %>% 
    group_by(a) %>% summarize(count=n())
  
  readl = cbind(rep(nt, nrow(readl)), rep(age, nrow(readl)), readl)
  names(readl) = c('name', 'age', 'readlength', 'count')
  

  coll.readlength = rbind(coll.readlength, readl)
  
  coll.ctp = rbind(coll.ctp, ctp)
}
options(scipen=999, digits=1) #shut off scientific notation

(plot_readlength = 
    ggplot(data=coll.readlength) +
    geom_col(aes(x=as.numeric(as.character(readlength)), y=count),
             colour='grey30') +
    facet_grid(~name) +
    xlab("Read Length (bp)") + 
    ylab("Count") + 
    theme_linedraw() +
    theme(
      axis.text = element_text(size = 8),
      axis.text.x  = element_text(angle = 70, hjust=1),
      axis.text.y  = element_text(angle = 70, hjust=1),
      axis.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      strip.text = element_text(size=6),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank()
    )
    )


(plot_plast = ggplot(data=coll.ctp) +
  geom_smooth( aes(x=z, y=CT), 
               lwd=0.5,
               alpha=0.8, 
               se=T,
               colour='black', 
               method = "gam", formula = y ~ s(x, k = 4)) +
  geom_line( aes(x=z, y=CT), alpha=0.5, colour='steelblue') +
  facet_grid(~name) + 
    xlim(0,100) + 
    ylim(0,0.10) + 
    xlab("Distance from 5' end (bp)") +
    ylab("C to T Frequency") +
    theme_linedraw() + 
    theme(
    axis.text = element_text(size = 8),
    axis.text.x  = element_text(angle = 70, hjust=1),
    axis.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    strip.text = element_text(size=6),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )
  )

gsa = ggarrange(plot_plast, plot_readlength, ncol=1, nrow=2, labels="AUTO")
ggsave(filename='damage_fig.png', 
       plot=gsa, 
       device=NULL, 
       width = 7.25, 
       height=3.75, 
       dpi=600)




  