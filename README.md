# LMTAG pipeline: GWAS leveraging local genetic correlation

## Step1: Run SUPERGNOVA

Default traits: 

"BC","LC-Adeno","LC-All","LC-Squam","PC","Asthma","Crohn","Eczema","IBD","PBC","TC","RA","SLE","UC","AD" ,"ADHD"  ,"ADs","AN" ,"ASD","BD","CP","CUD","DrnkWk","TG","epilepsy","Insomnia" ,"MDD","NSM","OCD","T2D","SA","SCZ","SD","SmkInit","TS","2HG","BMI","CAD","CKD","FG_adj","FI_adj","FP" ,"HBA1C" ,"HDL","Height","HOMA-B","HOMA-IR"  ,"HR" ,"IS" ,"LDL"  

The format of input trait GWAS sumstats should meet the requirements of munge function in LDSC!

1. Run SUPERGNOVA between above traits and primary trait. Commend for running SUPERGNOVA:

```
sh /ysm-gpfs/pi/zhao/yz738/BGM/pipeline.sh 54000 541000 AD BD  file1 file2  output_folder_name
```

2. create pipeline commend files. It will automatically generate folders and commend scripts for running SUPERGNOVA. ex:

```
sh /ysm-gpfs/pi/zhao/yz738/BGM/pipeline/pipeline.sh 66756 228951 IPF_IND_Finngen BC /ysm-gpfs/pi/zhao-data/yz738/GWAS/LC/Adeno/Adeno.txt /ysm-gpfs/pi/zhao/yz738/GWAS/BC/BC.txt /ysm-gpfs/pi/zhao/mc2792/supergnova/IPF_IND_Finngen_BC
``` 

3. create step commend files to run job files generated for pipeline.sh in batch. The order of commend is: run_munge.sh, run_step1.sh,run_step2.sh,run_step3.sh,run_step4.sh. example:

```
cd /ysm-gpfs/pi/zhao/mc2792/supergnova/IPF_IND_Finngen_BC/dsq
sh run_munge.sh
cd /ysm-gpfs/pi/zhao/mc2792/supergnova/IPF_IND_Finngen_BC/dsq
sh run_step1.sh
...
...
...
```

4. Organize the output format to the input of LMTAG using /ysm-gpfs/pi/zhao/mc2792/supergnova/supergnova_cmd/move.sh

5. FDR control on significant regions using ```ipf_allregion.R```

## Step 2: Run LMTAG

1. create region file(from ipf_allregion.R). create folder to storge the output of SUPERGNOVA on selected traits 

2. run LMTAG pipline (/ysm-gpfs/pi/zhao/mc2792/lmtag/code)

Format: sh pipeline.sh trait region data output. example: 

```
sh pipeline.sh IPF /ysm-gpfs/pi/zhao/mc2792/lmtag/example/ipf_region_sig.txt /ysm-gpfs/pi/zhao/mc2792/lmtag/example/data /ysm-gpfs/pi/zhao/mc2792/lmtag/output2
```

## Step 3: Post GWAS analysis


## Simulations





