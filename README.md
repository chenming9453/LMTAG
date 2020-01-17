# LMTAG
GWAS leveraging local genetic correlation

LMTAG pipline:

Step1: Run SUPERGNOVA

Default traits: 
"BC","LC-Adeno","LC-All","LC-Squam","PC","Asthma","Crohn","Eczema","IBD","PBC","TC","RA","SLE","UC","AD" ,"ADHD"  ,"ADs","AN" ,"ASD","BD","CP","CUD","DrnkWk","TG","epilepsy","Insomnia" ,"MDD","NSM","OCD","T2D","SA","SCZ","SD","SmkInit","TS","2HG","BMI","CAD","CKD","FG_adj","FI_adj","FP" ,"HBA1C" ,"HDL","Height","HOMA-B","HOMA-IR"  ,"HR" ,"IS" ,"LDL"  

1. Run SUPERGNOVA between above traits and primary trait
2. FDR control on significant regions using ipf_allregion.R 

Step 2: Run LMTAG
Format: sh pipeline.sh trait region data output
example: sh pipeline.sh IPF /ysm-gpfs/pi/zhao/mc2792/lmtag/example/ipf_region_sig.txt /ysm-gpfs/pi/zhao/mc2792/lmtag/example/data /ysm-gpfs/pi/zhao/mc2792/lmtag/output2

