#modified by Ming 09/30/2019
library(data.table)
options(stringsAsFactors=F)
args = commandArgs(trailingOnly=TRUE)

chr = as.numeric(args[1]) #chromosome number
munge1 <- as.character(args[2]) #gwas summary file
Dir <- as.character(args[3]) #output file 

ref_bim = paste0('/ysm-gpfs/pi/zhao/yz738/1000G/eur_chr', chr, '_SNPmaf5.bim')
file1_bim = read.table(gzfile(munge1), header=T) #read in gwas summary 
data.1kg = read.table(gzfile(ref_bim), header=F)
data.1kg = data.1kg[! data.1kg[,2] %in% data.1kg[duplicated(data.1kg[,2]),2],]
overlap = intersect(file1_bim[,1], data.1kg[,2]) #overlap between ref_bim and gwas summary 
data.1kg = data.1kg[data.1kg[,2] %in% overlap,]
data.1kg = data.1kg[! ((data.1kg[,5]=="G")&(data.1kg[,6]=="C") | (data.1kg[,5]=="C")&(data.1kg[,6]=="G") | (data.1kg[,5]=="A")&(data.1kg[,6]=="T") | (data.1kg[,5]=="T")&(data.1kg[,6]=="A")),]
output = data.frame(data.1kg[,c(1,4,2,3)])

#write _ovp.txt file 
output.whole = data.frame(data.1kg[,c(1,4,2,3)], rep(1,nrow(data.1kg)))
write.table(output[,3], paste0(Dir, '/step1/Ref/chr', chr, '_ovp.txt'), quote=F, row.names=F, col.names=F)

#add chr_block columns to output
chr_block = paste0('/ysm-gpfs/pi/zhao/yz738/Partition/block/chr', chr, '.txt')
block.file = read.table(chr_block, head=T)
temp <- c()
annotation <- paste0('BLOCK', 1:nrow(block.file))
for(j in 1:nrow(block.file)){
  start.snp <- block.file$start[j]
  stop.snp <- block.file$stop[j]
  temp.col <- output[,2]>=start.snp & output[,2]<=stop.snp
  temp <- cbind(temp, temp.col)
}
output <- cbind(output, temp)
gz = gzfile(paste0(Dir, '/step1/annotation/chr', chr,'.annot.gz'), "w") #annotation for region
write.table(output, gz, quote=F, col.names=c('CHR','BP','SNP','CM', annotation), row.names=F)
close(gz)
