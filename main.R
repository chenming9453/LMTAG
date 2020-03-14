
library(data.table)
options(stringsAsFactors=F)
args = commandArgs(trailingOnly=TRUE)

trait <- as.character(args[1]) #given major trait name
region <- as.character(args[2]) #region file
data_path <- as.character(args[3]) #given data path,which contains folders from SUPERGNOVA results
output <- as.character(args[4])

#Read in snplist 
snplist = fread("/ysm-gpfs/pi/zhao/mc2792/lmtag/code/snplist.txt")
traits_list = unlist(strsplit(list.files(data_path),"_"))
traits = unique(traits_list) #trait list
T = length(traits)
trait_o = traits[traits != trait]

#file path for each sumstats
get_file_path_sumstat = function(data_path){
   file_path_sumstat = data.frame(traits = traits,path = rep(NA,length(traits)))
   for(i in list.files(data_path)){
      name = strsplit(i,"_")
      path1 = paste0(data_path,"/",i,"/Munge/",name[[1]][1],".txt.sumstats.gz")
      file_path_sumstat$path[which(file_path_sumstat$traits == name[[1]][1])] = path1
      if(all(!is.na(file_path_sumstat$path))){
         break
      }
      path2 = paste0(data_path,"/",i,"/Munge/",name[[1]][2],".txt.sumstats.gz")
      file_path_sumstat$path[which(file_path_sumstat$traits == name[[1]][2])] = path2
      if(all(!is.na(file_path_sumstat$path))){
         break
      }
   }
   return(file_path_sumstat)
}
file_path_sumstat = get_file_path_sumstat(data_path)

step1.joblist <- c()
step2.joblist <- c()
for(s in 1:nrow(file_path_sumstat)){
   trait = file_path_sumstat$traits[s]
   munge = file_path_sumstat$path[s]
   Dir = paste0(output,"/",trait)
   #Step1
   for(i in 1:22){
     #use step1.R to generate /step1/annotation/chri.annot.gz and /step1/Ref/chri_ovp.txt
     annotation.joblist <- c(paste0('if [ ! -f "', Dir, '/step1/annotation/chr', i, '.annot.gz"', ' ]; then'),
                             paste0('Rscript --vanilla /ysm-gpfs/pi/zhao/mc2792/lmtag/code/step1.1.R ', i, ' ', munge, ' ', Dir),
                             paste0('fi'))
     #create annotation to /step1/Ref/chri
     annotation.joblist <- c(annotation.joblist,
                             paste0('/ysm-gpfs/pi/zhao/yz738/software/plink --bfile /ysm-gpfs/pi/zhao/yz738/1000G/eur_chr', i, '_SNPmaf5 --extract ', Dir, '/step1/Ref/chr', i, '_ovp.txt --make-bed --out ', Dir, '/step1/Ref/chr', i),

                             #calculate LDscore as /LDscore/ldsc/LDscore.i.l2.ldscore.gz
                             paste0('if [ ! -f "', Dir, '/LDscore/ldsc/LDscore.', i, '.l2.ldscore.gz"', ' ]; then'),
                             paste0('python2 /ysm-gpfs/pi/zhao/yz738/software/ldsc/ldsc/ldsc.py --l2 --bfile ', Dir, '/step1/Ref/chr', i, ' --ld-wind-cm 1 --annot ', Dir, '/step1/annotation/chr', i, '.annot.gz --out ', Dir, '/LDscore/ldsc/LDscore.', i),
                             paste0('fi'),

                             #re-format LDscore as /LDscore/region/chri.blockldscore.gz
                             paste0('if [ ! -f "', Dir, '/LDscore/region/chr', i, '.blockldscore.gz"', ' ]; then'),
                             paste0('Rscript --vanilla /ysm-gpfs/pi/zhao/mc2792/lmtag/code/step1.2.R ', i, ' ', Dir),
                             paste0('fi'))

     write.table(annotation.joblist, paste0(Dir, '/code/step1/step1chr', i, '.sh'), quote=F, row.names=F, col.names = F)
     step1.joblist <- c(step1.joblist, paste0('sh ', Dir, '/code/step1/step1chr', i, '.sh'))
   }
   #Step2
   for(i in 1:22){
     step2.temp <- c(paste0('if [ ! -f "', Dir, '/step2/chr', i, '.txt"', ' ]; then'),
                     paste0('Rscript --vanilla /ysm-gpfs/pi/zhao/mc2792/lmtag/code/step2.R ', i, ' ', munge, ' ', Dir),
                     paste0('fi'))
     write.table(step2.temp, paste0(Dir, '/code/step2/step2chr', i, '.sh'), quote = F, row.names = F, col.names = F)
     step2.joblist <- c(step2.joblist, paste0('sh ', Dir, '/code/step2/step2chr', i, '.sh'))
   }
   
}#for each trait
write.table(step1.joblist, paste0(output, '/code/step1.sh'), quote=F, row.names = F, col.names = F)
run.step1 <- c('source ~/.bashrc;', 'module load GCC;', 'partition=scavenge;', 'module load dSQ', paste0('dSQ --jobfile ', output, '/code/step1.sh -p ${partition} -n 1 -C avx2 --mem-per-cpu=12g -t 24:00:00 --mail-type=ALL --batch-file step1.pbs'), 'sbatch step1.pbs')
write.table(run.step1, paste0(output, '/code/run_step1.sh'), quote=F, row.names=F, col.names = F)

write.table(step2.joblist, paste0(output, '/code/step2.sh'), quote = F, row.names = F, col.names = F)
run.step2 <- c('source ~/.bashrc;', 'module load GCC;', 'partition=scavenge;', 'module load dSQ', paste0('dSQ --jobfile ', output, '/code/step2.sh -p ${partition} -n 1 -C avx2 --mem-per-cpu=10g -t 24:00:00 --mail-type=ALL --batch-file step2.pbs'), 'sbatch step2.pbs')
write.table(run.step2, paste0(output, '/code/run_step2.sh'), quote = F, row.names = F, col.names = F)



