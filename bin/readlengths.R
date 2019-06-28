lis = list.files('~/NeotomaSeq/results', recursive=T, pattern='rmdup.fq', full.names = T)

full_table = read.delim("./Files2Run/packrat_aDNA.csv", sep = ",",  na.strings = "null")

nam = list.dirs('~/NeotomaSeq/results/')[2:12]
readlengths <- function(x){
  targ = strsplit(nam[x], '/')[[1]][7]
  #comm = paste("awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' ", lis[x], " > ", targ, ".readlengths", sep = '')
  comm = paste("awk 'NR%4 == 2 {print length($0)}' ", lis[x], " > ", targ, ".readlengths", sep = '')
  print(comm)
  system(comm)
}

#readlengths(1)

n_cores <- 11
cl <- makeCluster(n_cores)
clusterExport(cl=cl, varlist=c("nam", "lis"))
parSapply(cl,1:11,FUN=readlengths)
stopCluster(cl)

readlengths = list.files(pattern="*.readlengths")
df = data.frame(length=numeric(),  extraction=character())
for(j in 1:length(readlengths)){
  sub = data.table::fread(readlengths[j])
  sub$V3 = rep(readlengths[j], nrow(sub))
  df = base::rbind(df, sub, fill=T)
}

df = df[,3:4]
colnames(df) = c('length', 'sample')

df = df[sample(1:nrow(df), 50000),]

#Correct sample names
df$Age = rep(NA, nrow(df))

for(zz in 1:nrow(df)){ 
  df$sample[zz] = strsplit(df$sample[zz], '[.]')[[1]][1]
  name = strsplit(df$sample[zz], '[-]')[[1]][1]
  if(name == "GC100B"){next} #Make another figure to compare MX and US sites
  df$Age[zz] = as.numeric(as.character(unique(full_table[grep(name, full_table$Name), 'Midden.age'])))
 
   

}

df = na.omit(df)
#sample name as factor ordered by age

name.map=unique(df[,c('sample', 'Age')])
uni = unique(name.map[order(name.map$Age),]$sample)
df$sample = factor(x=df$sample, ordered=T, levels=uni)



rl1 = ggplot(df) + 
  geom_density(aes(df$length)) + 
  facet_grid(~sample, scales = 'free_x', space='free') + 
  theme_linedraw() +
  theme(    axis.text = element_text(size = 12),
            axis.text.x  = element_text(angle = 90, size =12),
            axis.title = element_text(size = 12),
            legend.text = element_text(size = 12)
            ) +
  xlab('Merged Read Length (bp)') + 
  ylab('Density Estimation') + xlim(c(20, 250))
rl1


png(filename = 'readlengths.png', height = 4, width = 10.25, units = 'in', res=1200)
plot_grid(rl1)
dev.off()
