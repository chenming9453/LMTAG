#modified by Ming 09/26/2019
library(data.table)
options(stringsAsFactors=F)
args = commandArgs(trailingOnly=TRUE)

i = as.numeric(args[1]) #chr
Dir = as.character(args[2]) #output path

annot.file <- paste0(Dir, '/step1/annotation/chr', i,'.annot.gz') #output from step 1.1
annot <-  read.table(gzfile(annot.file), head=T)

ldscore.file <- paste0(Dir, '/LDscore/ldsc/LDscore.', i, '.l2.ldscore.gz') #from SUPERGNOVA.R output from LDSC format
ldscore <- read.table(gzfile(ldscore.file), head=T)

nr <- nrow(annot)
output <- data.frame('CHR'=numeric(), 'SNP'=character(), 'BP'=numeric(), 'L2'=numeric())
nc.annot <- ncol(annot)
nc.ldscore <- ncol(ldscore)
pb = txtProgressBar(0, nr, style=3)
for(j in 1:nr){
  annot.sub <- as.numeric(annot[j, 5:nc.annot])
  ldscore.sub <- as.numeric(ldscore[j, 4:nc.ldscore])

  l2 <- t(annot.sub)%*%ldscore.sub
  temp <- data.frame('CHR'=i, 'SNP'=ldscore[j,2], 'BP'=ldscore[j,3], 'L2'=l2)
  output <- rbind(output, temp)
  setTxtProgressBar(pb, j)
}
gz = gzfile(paste0(Dir, '/LDscore/region/chr',i,'.blockldscore.gz'), 'w') #ouput region LDSC 
write.table(output, gz, quote=F, row.names=F)
close(gz)
