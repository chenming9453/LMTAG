#modified by Ming 09/30/2019 calculating heritability
library(data.table)
options(stringsAsFactors=F)
args = commandArgs(trailingOnly=TRUE)
chr = as.numeric(args[1]) #chromosome  
munge = as.character(args[2]) #summary statistics AD.txt.sumstats.gz
Dir =as.character(args[3]) #output directory /AD

snp_region = fread("/ysm-gpfs/pi/zhao/mc2792/lmtag/code/snplist.txt")
sumstats1 = fread(munge)
N1 = sumstats1$N[1]
LDscore.region <- fread(paste0(Dir, '/LDscore/region/chr', chr, '.blockldscore.gz')) #reformated ld score


overlap.SNP = intersect(sumstats1$SNP,LDscore.region$SNP)  #370874 for chr 1
snp_region_sub = snp_region[snp_region$snp %in% overlap.SNP]

snp_region_sub = snp_region_sub[snp_region_sub$chr == chr,]
snp_region_list = snp_region_sub[, ID := .GRP, by = .(start,end)]
region_index = unique(snp_region_sub[,c('start','end')])

region_num = nrow(region_index) #2353
h = rep(0,region_num)
m = rep(0,region_num)
for(i in 1:region_num){
    snp = snp_region_list$snp[snp_region_list$ID == i] #2913
    m[i] <- length(snp)  
    sum1_sub <- sumstats1[sumstats1$SNP %in% snp,] #2337
    LDscore.block <- LDscore.region[LDscore.region$SNP %in% snp,] #2337
    m[i] = nrow(sum1_sub)
    L2 = as.numeric(as.character(LDscore.block$L2))
    h[i] = m[i]* (mean((sum1_sub$Z)^2)-1)/(N1*mean(L2))
}
region_index$h = h
region_index$m = m
write.table(region_index, file=paste0(Dir, '/step2/chr', chr, '.txt'), quote = FALSE, sep ="\t", row.names = F)

