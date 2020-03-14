library(dplyr)
library(data.table)
library(MASS)
library(mixtools)
library(fitdistrplus)
options(stringsAsFactors=F)
args = commandArgs(trailingOnly=TRUE)
trait_s <- as.character(args[1]) #given major trait name
region_file <- as.character(args[2]) #region file
data_path <- as.character(args[3]) #given data path,which contains folders from SUPERGNOVA results
output <- as.character(args[4])
snplist = fread("/ysm-gpfs/pi/zhao/mc2792/lmtag/code/snplist.txt")


trait_s = "IPF"
region_file = "/ysm-gpfs/pi/zhao/mc2792/supergnova_cmd/sig_region_allen_new.txt"
data_path = "/ysm-gpfs/pi/zhao/mc2792/lmtag/IPF-Allen-NEW"
output = "/ysm-gpfs/pi/zhao/mc2792/lmtag/IPF-Allen-NEW-OUTPUT"


trait_s = "IPF-IND-Fingerlin"
region_file = "/ysm-gpfs/pi/zhao/mc2792/supergnova_cmd/sig_region_fingerlin_new"
data_path = "/ysm-gpfs/pi/zhao/mc2792/lmtag/IPF-Fingerlin-NEW-2"
output = "/ysm-gpfs/pi/zhao/mc2792/lmtag/IPF-Fingerlin-NEW-2-OUTPUT"
traits = unique(unlist(strsplit(list.files(data_path),"_"))) #trait list: 37 traits + IPF
T = length(traits)

#####################################################################################################################
#1. read in significant region
region = read.table(region_file,header = T)[,1:4] #87 significant regions
region_unique = distinct(region, start,chr, .keep_all= TRUE) #50 unique regions
region_unique$index = seq(1,nrow(region_unique),1)
region_unique  = subset(region_unique,select = c(index, chr,start,end))
index_trait = list()
for(i in 1:nrow(region_unique)){
  index_trait[[i]] =subset(region,region$chr == region_unique$chr[i] & region$start == region_unique$start[i])$trait
}

#####################################################################################################################
#2. read in intercept fo LDSC to create Sigma_LMTAG: 
get_region_sigma_lmtag = function(region_traits){
  temp = matrix(0,nrow = length(region_traits),ncol = length(region_traits) ,dimnames = list(region_traits,region_traits))
  #diag
  for(i in 1:length(region_traits)){
    trait = region_traits[i]
    tempdata = fread(file_path_sumstat$path[file_path_sumstat$traits == trait])
    N = tempdata$N[1]
    temp[i,i] = 1/N
  }
  #off_diag
  for(i in 1:(length(region_traits)-1)){
    for(j in (i+1) :length(region_traits)){
      t1 = region_traits[i]
      t2 =region_traits[j]
      #read intercept
      suppressWarnings({
        filename = paste0(data_path,"/",t1,"_",t2,"/Whole/",t1,"_",t2,".txt")
        file = try(fread(filename), silent=TRUE)
        if(class(file) == "try-error"){
          filename = paste0(data_path,"/",t2,"_",t1,"/Whole/",t2,"_",t1,".txt")
        }
        inter = fread(filename)
      })#suppressWarnings
      temp[i,j] = inter$intercept
      temp[j,i] = temp[i,j]
    }
    
  }
  return(temp)
}

###########################################################################################
#3. read in local genetic covariance between trait s and other traits
get_local_rho = function(region){
  local_rho = data.frame(chr = numeric(),start = numeric(),end = numeric(),trait = character(),rho = numeric())
  suppressWarnings({
    for(i in 1:nrow(region)){
      file_name = paste0(data_path,"/",trait_s,"_",region$trait[i],"/Local/chr",region$chr[i],".txt")
      file = try(fread(file_name), silent=TRUE)
      if(class(file) == "try-error"){
        file_name = paste0(data_path,"/",region$trait[i],"_",trait_s,"/Local/chr",region$chr[i],".txt")
      }
      data = fread(file_name)
      temp_region = paste0("tilde.chr",region$chr[i],".",region$start[i],".",region$end[i],".gz.txt")
      temp_data = data.frame(chr = region$chr[i],start = region$start[i],end = region$end[i], trait = region$trait[i], rho = data$rho[data$region == temp_region])
      local_rho  = rbind(local_rho,temp_data)
    }
  })
  local_rho = distinct(local_rho)
  return(local_rho)
}
local_rho = get_local_rho(region)

###########################################################################################
#4. read in gwas beta for t traits, only obtain the snps for sig regions
#####get the file path for each trait's sumstat 
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

###get the beta for each snp 
get_snp_beta = function(region){
  snp_beta = data.frame(chr = numeric(),start = numeric(),end = numeric(), snp = character(),beta = character(),trait = numeric())
  for(i in unique(region$trait)){
    #subset snplist
    region_sub = subset(region,region$trait==i)
    sub_list = data.frame()
    for(j in 1:nrow(region_sub)){
      temp = subset(snplist, chr == region_sub$chr[j] & start >=  region_sub$start[j] & end <= region_sub$end[j])
      sub_list = rbind(sub_list,temp)
    }
    tempdata = fread(file_path_sumstat$path[file_path_sumstat$traits == i])
    data_merge = merge(tempdata,sub_list,by.x = "SNP",by.y = "snp")
    data_merge$BETA  = data_merge$Z/sqrt(data_merge$N)
    data_merge = subset(data_merge, select = c("chr","start","end","SNP","BETA"))
    data_merge$trait = rep(i,nrow(data_merge))
    colnames(data_merge) = c("chr","start","end","snp","beta","trait")
    snp_beta = rbind(snp_beta,data_merge)
  }
  sub_list = data.frame()
  for(i in 1:nrow(region_unique)){
    temp = subset(snplist, chr == region_unique$chr[i] & start >=  region_unique$start[i] & end <= region_unique$end[i])
    sub_list = rbind(sub_list,temp)
  }
  #for trait_s
  tempdata = fread(file_path_sumstat$path[file_path_sumstat$traits == trait_s])
  data_merge = merge(tempdata,sub_list,by.x = "SNP",by.y = "snp")
  data_merge$BETA  = data_merge$Z/sqrt(data_merge$N)
  data_merge = subset(data_merge, select = c("chr","start","end","SNP","BETA"))
  data_merge$trait = rep(trait_s,nrow(data_merge))
  colnames(data_merge) = c("chr","start","end","snp","beta","trait")
  snp_beta = rbind(snp_beta,data_merge)
  return(snp_beta)
}
snp_beta = get_snp_beta(region)

#####################################################################################################################
#5. for Omeg_j off diagonal elements
get_omeg_off = function(region_unique,snplist){
  snp_l2 = list()
  for(i in 1:nrow(region_unique)){
    tempdata = data.frame(snp = character(),trait1 = character(),trait2 = character(),L2 = numeric(),rho_m = numeric())
    temp_rho_m = data.frame(trait1 = character(),trait2 = character(),rho_m = numeric())
    chr = region_unique$chr[i]
    start = region_unique$start[i]
    end = region_unique$end[i]
    temp = which(snplist$chr == chr & snplist$start == start & snplist$end == end)
    subsnp = snplist[temp,]
    trait = index_trait[[i]]
    trait = c(trait,trait_s)
    for(j in 1:(length(trait)-1)){
      for(k in (j+1):length(trait)){
        trait1 = trait[j]
        trait2 = trait[k]
        ####get ldscore for each snp
        suppressWarnings({
          filename = paste0(data_path,"/",trait1,"_",trait2,"/LDscore/region/chr",chr,".blockldscore.gz")
          file = try(fread(filename), silent=TRUE)
          if(class(file) == "try-error"){
            filename = paste0(data_path,"/",trait2,"_",trait1,"/LDscore/region/chr",chr,".blockldscore.gz")
          }
          temp_l2 = fread(filename)
        })#suppressWarnings
        temp_l2_snp = merge(subsnp,temp_l2,by.x = "snp",by.y = "SNP")  
        temp_l2_snp$trait1 = trait1
        temp_l2_snp$trait2 = trait2
        temp_l2_snp = subset(temp_l2_snp,select = c("snp","trait1","trait2","L2")) 
        
        ####get rho over m 
        suppressWarnings({
          filename_2 = paste0(data_path,"/",trait1,"_",trait2,"/Local/chr",chr,".txt")
          file_2 = try(fread(filename_2), silent=TRUE)
          if(class(file_2) == "try-error"){
            filename_2 = paste0(data_path,"/",trait2,"_",trait1,"/Local/chr",chr,".txt")
          }
          temp = fread(filename_2)
        })#suppressWarnings
        region_name = paste0("tilde.chr",chr,".",start,".",end,".gz.txt")
        rho_m_temp = temp$rho[which(temp$region == region_name)]/temp$m1[which(temp$region == region_name)]
        temp_l2_snp$rho_m = rho_m_temp
        
        tempdata = rbind(tempdata,temp_l2_snp)
      }#trait s and trait t
    }
    snp_l2[[i]] = tempdata
  }
  return(snp_l2)
}
omeg_off = get_omeg_off(region_unique,snplist)

########################################################################################################
#6. get regional h2, l2, m for each trait 
###first need to generate the files according to SUPERGNOVA.R step1.1 step1.2 and step3
get_omeg_diag = function(region_unique,snplist){
  h_m = data.frame(chr = numeric(),start = numeric(),end = numeric(), trait = character(),h = numeric(),m = numeric())
  snp = data.frame(chr = numeric(),start = numeric(),end = numeric(), trait = character(),snp = character(),l2 = numeric())
  for(i in 1:nrow(region_unique)){
    temp = which(snplist$chr == region_unique$chr[i] & snplist$start == region_unique$start[i] & snplist$end == region_unique$end[i])
    sub_snp = snplist[temp,]$snp
    traits = index_trait[[i]]
    traits = c(trait_s,traits)
    chr = region_unique$chr[i]
    start = region_unique$start[i]
    end = region_unique$end[i]
    temp_h = c()
    temp_m = c()
    for(j in traits){
      #update h_m
      file = paste0(output,"/",j,"/step2/chr",chr,".txt")
      temp_data = fread(file)
      temp_h = c(temp_h,temp_data$h[which(temp_data$start == start & temp_data$end == end)])
      temp_m = c(temp_m,temp_data$m[which(temp_data$start == start & temp_data$end == end)])
      #update snp
      file_l2 = paste0(output,"/",j,"/LDscore/region/chr",chr,".blockldscore.gz")
      temp_l2 = fread(file_l2)
      temp_l2_snp = temp_l2[temp_l2$SNP %in% sub_snp]
      temp_l2_snp$chr = chr
      temp_l2_snp$start = start
      temp_l2_snp$end = end 
      temp_l2_snp$trait = j
      temp_snp = subset(temp_l2_snp,select = c("chr","start","end","trait","SNP","L2"))
      colnames(temp_snp) = c("chr","start","end","trait","snp","l2")
      snp = rbind(snp,temp_snp)
    }
    temp_h_m = data.frame(chr = chr, start = start, end = end ,trait = traits, h = temp_h, m = temp_m)
    h_m = rbind(h_m,temp_h_m)
  }
  omega_diag = list("h_m" = h_m, "snp" = snp)
  return(omega_diag)
}
omeg_diag = get_omeg_diag(region_unique,snplist)

########################################################################################################
# 7.construct Omega_j    
get_omeg = function(omeg_off,omeg_diag,region,snp){
  region_traits = index_trait[[region]]
  region_traits = c(trait_s,region_traits)
  N = length(region_traits)
  chr = region_unique$chr[region]
  start = region_unique$start[region]
  end = region_unique$end[region]
  h_m = omeg_diag$h_m
  off = omeg_off[[region]]
  #regional l2 ffor snps
  snp_l2 = omeg_diag$snp
  temp = which(snp_l2$chr == chr & snp_l2$start == start)
  snp_l2_region = snp_l2[temp,]
  
  omega_j = matrix(0,nrow = N,ncol = N,dimnames = list(region_traits,region_traits))
  #diag of omega_j
  for(i in 1:N){
    trait = region_traits[i]
    index = which(h_m$trait == trait & h_m$chr == chr & h_m$start == start)
    omega_j[i,i] = snp_l2_region$l2[which(snp_l2_region$snp == snp & snp_l2_region$trait == trait)] * h_m$h[index]/h_m$m[index]
  }
  #off-diag of omega_j
  for(i in 1:(N-1)){
    for(j in (i+1):N){
      t1 = region_traits[i]
      t2 = region_traits[j]
      index = which(off$snp ==snp & off$trait1 == t1 & off$trait2 == t2)
      if(length(index) == 0){
        index = which(off$snp==snp & off$trait1 == t2 & off$trait2 == t1)
      }
      omega_j[i,j] = off[index,]$L2 * off[index,]$rho_m
      omega_j[j,i] = omega_j[i,j]
    }
  }
  return(omega_j)
}

########################################################################################################
#8.test
##for each region:
##calculate W_j,s 
##calculate beta/se/z/p/
get_beta = function(omeg_diag,threshold){
  beta = data.frame(index = numeric(),chr = numeric(),start = numeric(),end = numeric(),snp = character(),beta = numeric(),var = numeric(), se = numeric(), z = numeric(),p = numeric())
  omeg_diag_h_m = omeg_diag[['h_m']]
  omeg_diag_snp = omeg_diag[['snp']]
  
  eigen_data = data.frame(snp = character(),beta = numeric(),var = numeric(), min_eigen=numeric(),max_eigen = numeric())
  
  for(i in 1:nrow(region_unique)){
    chr = region_unique$chr[i]
    start = region_unique$start[i]
    end = region_unique$end[i]
    region_traits = c(trait_s,index_trait[[i]])
    rsig_lmtag = get_region_sigma_lmtag(region_traits)
    h2_rjs = omeg_diag_h_m$h[which(omeg_diag_h_m$trait == trait_s & omeg_diag_h_m$chr == chr & omeg_diag_h_m$start == start )]
    m_rjs = omeg_diag_h_m$m[which(omeg_diag_h_m$trait == trait_s & omeg_diag_h_m$chr == chr & omeg_diag_h_m$start == start )]
    #order rho_rjs
    rho_rjs = rep(h2_rjs,length(region_traits))
    names(rho_rjs) = region_traits
    for(s in index_trait[[i]]){
      index = which(local_rho$chr == chr & local_rho$start == start & local_rho$trait == s)
      rho_rjs[which(names(rho_rjs)==s)] = local_rho$rho[index]
    }
    rho_rjs[1] = h2_rjs
    rho_h2 = rho_rjs/h2_rjs
    
    #get the overlapped local snp list ampng traits
    temp_snp_list = omeg_off[[i]]
    overlap_snp = temp_snp_list$snp
    if(length(region_traits)> 2){
      for(k in 1:(length(region_traits)-1)){
        for(t in (k+1):length(region_traits)){
          T1 = region_traits[k]
          T2 = region_traits[t]
          index = which(temp_snp_list$trait1==T1 & temp_snp_list$trait2 == T2)
          if(length(index) == 0){index = which(temp_snp_list$trait1==T2 & temp_snp_list$trait2 == T1)}
          temp_snp = temp_snp_list$snp[index]
          overlap_snp = intersect(overlap_snp,temp_snp)
        }
      }
      local_snp_list = overlap_snp 
    }else{local_snp_list = temp_snp_list$snp}
    
    for(j in local_snp_list){
      ####get beta_gwas_j#####
      beta_gwas_j = rep(1,length(region_traits))
      names(beta_gwas_j) = region_traits
      for(s in region_traits){
        beta_gwas_j[which(names(beta_gwas_j) == s)] = snp_beta[which(snp_beta$snp == j & snp_beta$trait == s)]$beta
      }
      beta_gwas_j = as.numeric(beta_gwas_j)
      
      #calculate A_js
      omeg_j = get_omeg(omeg_off,omeg_diag,i,j)  
      l_js = omeg_diag_snp$l2[which(omeg_diag_snp$snp == j & omeg_diag_snp$trait == trait_s)]
      rr_h2 = rho_rjs %*% t(rho_rjs)/ h2_rjs
      A_js = omeg_j + rsig_lmtag - l_js/m_rjs*rr_h2
      
      ###force small eigenvalue to be zero
      if(any(is.na(A_js))){
        var_j = NA
      }
      else{
        e = eigen(A_js)
        eigvalue = e$values
        if (any(eigvalue <= threshold)){
          var_j = NA
        }else{
          eigvalue = 1/eigvalue
          eigvalue[eigvalue == Inf] = 0
          D = diag(eigvalue)
          P = e$vectors
          out <- tryCatch(solve(P) %*% P, error = function(e) e)
          if(any(class(out) == "error")){
            W_js = NA
            beta_j = NA
            var_j = NA
          }else{
            W_js = P %*% D %*% solve(P)
            beta_j = t(rho_h2) %*% W_js %*% beta_gwas_j /(rho_h2 %*% W_js %*% rho_h2)
            var_j = 1/(rho_h2 %*% W_js %*% rho_h2)
            temp = data.frame(snp = j,beta = beta_j,var = var_j, min_eigen=min(e$values),max_eigen=max(e$values))
            eigen_data = rbind(eigen_data,temp)
            
          }
        }
      }#A_js contains NA
      
      if(is.na(var_j)){
        #print(paste0("chr",chr,"start",start,"end",end,"snp:",j))
      }else{
        if(var_j <= 0){temp_beta = data.frame(index = i, chr = chr, start = start, end = end, snp = j, beta = beta_j, var = var_j, se = "NA", z = "NA", p = "NA")
        }else{
          se_j = sqrt(var_j)
          z_j = beta_j/se_j
          p_j = 2 * pnorm(abs(z_j),lower.tail = F)
          temp_beta = data.frame(index = i, chr = chr, start = start, end = end, snp = j, beta = beta_j,var = var_j, se = se_j, z = z_j, p = p_j)
        }
        beta = rbind(beta,temp_beta)
      }#NA variance
    }#for snps in the region
    
    #all_egen[[i]] = temp
  }#region
  write.table(eigen_data,file = paste0("/ysm-gpfs/pi/zhao/mc2792/lmtag/eigen_new_",trait_s,"_",threshold,".txt"), sep =  "\t", row.names = FALSE,quote = FALSE)
  return(beta)
}

threshold = 0
beta = get_beta(omeg_diag,threshold)
write.table(beta,file = paste0("/ysm-gpfs/pi/zhao/mc2792/lmtag/lmtag_output_0_new_",trait_s,".txt"), sep =  "\t", row.names = FALSE,quote = FALSE)


