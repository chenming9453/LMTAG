import sys
import pandas as pd
from pandas.core.common import flatten
import numpy as np
import os
import scipy.stats
import math
import time
pd.set_option("display.max_colwidth",10000)
start = time.clock()


trait_1 = sys.argv[1]
trait_2 = sys.argv[2]
chr = int(sys.argv[3])
data_path = sys.argv[4]
output = sys.argv[5]


####test mode#####
data_path = "/home/mc2792/scratch60/lmtag/supergnova"
trait_1 = "ASD"
trait_2 = "ADHD"
chr = 1
output = "/home/mc2792/scratch60/lmtag/lmtag/output/"

snplist = pd.read_csv("/ysm-gpfs/pi/zhao/mc2792/lmtag/code/snplist.txt", delimiter = '\t')
traits = [trait_1,trait_2]
T = 2


#####Test mode######




#########################################################################################################################################################
########################################################################################################################################################################################
def get_file_path_sumstat(data_path):
    file_path_sumstat = pd.DataFrame({'traits':traits,'path':[None] * len(traits)})
    for i in arr:
        if not any(x is None for x in file_path_sumstat.path):
            break
        name = i.split("_")
        path1 = data_path + "/" + i + "/Munge/"+ name[0]+".txt.sumstats.gz"
        file_path_sumstat.loc[file_path_sumstat.traits == name[0],'path']=path1
        path2 = data_path + "/" + i + "/Munge/"+ name[1]+".txt.sumstats.gz"
        file_path_sumstat.loc[file_path_sumstat.traits == name[1], 'path'] = path2
    return file_path_sumstat

def beta_n(trait):
    temp = file_path_sumstat.loc[file_path_sumstat.traits == trait,'path']
    file_path = temp.to_string().split(" ")[4]
    gwas_1 = pd.read_csv(file_path, delimiter=' ')
    n1 = gwas_1.loc[0,'N']
    snp_gwas_1 = pd.merge(snplist,gwas_1,how = "inner", left_on='snp', right_on='SNP')
    snp_gwas_1['beta'] = snp_gwas_1['Z']/snp_gwas_1['N']**(1/2)
    snp_beta_1 = snp_gwas_1[['chr','start','end','snp','beta']]
    return n1,snp_beta_1

def get_intercept(trait_1,trait_2):
    filename = data_path + '/' + trait_1 + '_' + trait_2 + '/Whole/' + trait_1 + '_' + trait_2 + '.txt'
    if not os.path.exists(filename):
        filename = data_path + '/' + trait_2 + '_' + trait_1 + '/Whole/' + trait_2 + '_' + trait_1 + '.txt'
    inter = pd.read_csv(filename, delimiter='\t')
    return inter.loc[0,'intercept']

def get_ldsc(trait_1,trait_2,chr):
    filename = data_path + "/" + trait_1 + "_" + trait_2 + "/LDscore/region/chr" + str(chr) + ".blockldscore.gz"
    if not os.path.exists(filename):
        filename = data_path + "/" + trait_2 + "_" + trait_1 + "/LDscore/region/chr" + str(chr) + ".blockldscore.gz"
    ld = pd.read_csv(filename,delimiter=' ')
    ld.columns = ['chr','snp','bp','lj']
    return ld

def get_local_rho(region,chr,trait_1,trait_2):
    file_name = data_path + '/' + trait_1 + '_' + trait_2 + '/Local/chr' + str(chr) + '.txt'
    t1 = trait_1
    t2 = trait_2
    if not os.path.exists(file_name):
        file_name = data_path + '/' + trait_2 + '_' + trait_1 + '/Local/chr' + str(chr) + '.txt'
        t2 = trait_1
        t1 = trait_2
    rho = pd.read_csv(file_name,delimiter='\t')
    rho.dropna(inplace=True)
    rho['chr'] = chr
    rho[['tilde','c','start','end','gz','txt']] = rho.region.str.split(".",expand = True)
    local_rho = rho[['chr','start','end','rho','var.rho','h1','var1','h2','var2','m1']]
    if t1 == trait_1:
        local_rho.columns = ['chr','start','end','rho','rho_var','h2_t1','h2_var_t1','h2_t2','h2_var_t2','m']
    else:
        local_rho.columns = ['chr', 'start', 'end', 'rho', 'rho_var', 'h2_t2', 'h2_var_t2', 'h2_t1','h2_var_t1','m']
    local_rho['start'] = pd.to_numeric(local_rho['start'])
    local_rho['end'] = pd.to_numeric(local_rho['end'])
    return local_rho


####################################################################################################################################
#####################################################################################################################
#1. read in region file

region = pd.read_csv("/Users/mingchen/Documents/2020_spring/project_IPF_integrate_analysis/main/code/region_code.txt",delimiter = '\t')
trait_rest = traits[traits != trait_1]
print("read in region...")
file_path_sumstat = get_file_path_sumstat(data_path)

#First get the beta and n for traits
print("get beta...")
n1,snp_beta_1 = beta_n(trait_1)
n2,snp_beta_2 = beta_n(trait_2)

print("read in intercept...")
# intercept(n12*rho12/n1n2)
intercept = get_intercept(trait_1,trait_2)

print("read in local rho...")
#rho_R(j),1,h1_R(j),h2_R(j),var(rho_R(j)),var(h1_R(j)),var(h2_R(j)), m(R(j))
#notice some regions h2_IPF<0
local_rho = get_local_rho(region,chr,trait_1,trait_2)

print("read in ld...")
#ldsc for each snp
ld = get_ldsc(trait_1,trait_2,chr)


#for each snp update:
lmtag_snp = list(set(snp_beta_1['snp']) & set(snp_beta_2['snp']) & set(ld['snp']))
n_len = len(lmtag_snp)
beta = pd.DataFrame(columns=['index', 'chr', 'start','end','snp','beta','var','se','z','p'])
print("get beta of "+str(len(lmtag_snp))+" snps for "+trait_2+'...')
for i in range(n_len):
    if i % 10000 ==0:
        print(i)
    j = lmtag_snp[i]
    ld_temp = ld[ld.snp == j]
    ld_temp = ld_temp.reset_index()
    l_j = ld_temp.ix[0,'lj']
    temp_chr = ld_temp.ix[0,'chr']
    temp_bp = ld_temp.ix[0,'bp']
    region_temp = local_rho[(local_rho.chr == temp_chr) & (local_rho.start <= temp_bp) & (local_rho.end >= temp_bp)]
    region_temp = region_temp.reset_index()
    if not region_temp.empty:
        temp_chr = region_temp.ix[0,'chr']
        temp_start = region_temp.ix[0,'start']
        region_index = region[(region.chr == temp_chr) & (region.start == temp_start)].reset_index().ix[0,'index']
        h2_Rj_1 = region_temp.ix[0,'h2_t1']
        h2_Rj_2 = region_temp.ix[0,'h2_t2']
        rho_Rj_1_2 = region_temp.ix[0,'rho']
        rho_var_Rj_1_2 = region_temp.ix[0,'rho_var']
        h2_var_Rj_1 = region_temp.ix[0,'h2_var_t1']
        m_Rj = region_temp.ix[0,'m']

        # get W_j1
        var_beta = np.array([[1/n1, intercept], [intercept, 1/n2+l_j/m_Rj*max(0,h2_Rj_2-rho_Rj_1_2**2/h2_Rj_1)]])
        W_j1 = np.linalg.inv(var_beta)

        # get beta
        beta_j_1 = snp_beta_1[snp_beta_1.snp == j].reset_index().ix[0,'beta']
        beta_j_2 = snp_beta_2[snp_beta_2.snp == j].reset_index().ix[0,'beta']
        beta_gwas = np.array([beta_j_1, beta_j_2])
        rho_Rj_1 = np.array([h2_Rj_1,rho_Rj_1_2])
        rho_h = rho_Rj_1/h2_Rj_1
        beta_lmtag_j_1 = (rho_h @ W_j1 @ beta_gwas) / (rho_h @ W_j1 @ rho_h)

        # get variance
        eta = rho_Rj_1_2 / h2_Rj_1
        var_eta = rho_var_Rj_1_2 / (h2_Rj_1**2) + h2_var_Rj_1 * (rho_Rj_1_2**2) / (h2_Rj_1**4)
        var_lmtag_j_1 = 1/(rho_h @ W_j1 @ rho_h)+var_eta*(n1**2*n2**2*(beta_j_1-beta_j_2)**2)/((n1+n2*eta)**4)

        se = var_lmtag_j_1**(1/2)
        z = beta_lmtag_j_1 / se
        p = 2 * (1-scipy.stats.norm(0, 1).cdf(abs(z)))
        beta = beta.append({'index': region_index, 'chr':temp_chr, 'start':temp_start,'end': region_temp.ix[0,'end'],'snp':j,'beta':beta_lmtag_j_1,'var':var_lmtag_j_1,'se':se,'z':z,'p':p}, ignore_index=True)


print("start to write file...")
outname = output+'/LMTAG_'+trait_2+'.txt'
beta.to_csv(outname, sep="\t")
end = time.clock()
print("The total running time is:")
print(end-start)
