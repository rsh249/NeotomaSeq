##Build taxonomizr database
library(taxonomizr)
libdir='taxonomizr'
dir.create(libdir)
setwd(libdir)
getNamesAndNodes()
getAccession2taxid()
read.accession2taxid(list.files('.','accession2taxid$'),'accessionTaxa.sql')
print(paste('taxonomizr database built and located at', getwd(), sep=' '))
      