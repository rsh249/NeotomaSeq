##Build taxonomizr database
library(taxonomizr)
args=commandArgs(trailingOnly=TRUE)
targ=args[[1]];
libdir=targ
dir.create(libdir)
setwd(libdir)
getNamesAndNodes()
getAccession2taxid(types=c('nucl_gb'))
getAccession2taxid()
system("gunzip *.gz")
read.accession2taxid(list.files('.','accession2taxid'),'accessionTaxa.sql')
print(paste('taxonomizr database built and located at', getwd(), sep=' '))
      
