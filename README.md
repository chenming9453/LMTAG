# UTMOST commend:

mkdir sample_data/results_GTEx

UTMOST_path = /home/mc2792/UTMOST/

python2 joint_GBJ_test.py \
--weight_db $UTMOST_path/sample_data/weight_db_GTEx/ \
--output_dir $UTMOST_path/sample_data/results_GTEx/ \
--cov_dir $UTMOST_path/sample_data/covariance_joint/ \
--input_folder $UTMOST_path/sample_data/results/ \
--gene_info $UTMOST_path/intermediate/gene_info.txt \
--output_name ipf_4_joint \
--start_gene_index 1 \
--end_gene_index 17290

# LMTAG pipeline: GWAS leveraging local genetic correlation

## Step1: Run SUPERGNOVA

Default traits: 

"BC","LC-Adeno","LC-All","LC-Squam","PC","Asthma","Crohn","Eczema","IBD","PBC","TC","RA","SLE","UC","AD" ,"ADHD"  ,"ADs","AN" ,"ASD","BD","CP","CUD","DrnkWk","TG","epilepsy","Insomnia" ,"MDD","NSM","OCD","T2D","SA","SCZ","SD","SmkInit","TS","2HG","BMI","CAD","CKD","FG_adj","FI_adj","FP" ,"HBA1C" ,"HDL","Height","HOMA-B","HOMA-IR"  ,"HR" ,"IS" ,"LDL"  

Folders with SUPERGNOVA results:

```/ysm-gpfs/pi/zhao/yz738/SUPERGNOVA/Results```

```/home/mc2792/scratch60/supergnova2```

```/home/mc2792/scratch60/supergnova```

```/ysm-gpfs/pi/zhao/mc2792/supergnova/```

```/ysm-gpfs/pi/zhao-data/mc2792/supergnova```

```/home/mc2792/scratch60/supergnova```

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

5. The local genetic correlations estimated from SUPERGNOVA are located in ```/ysm-gpfs/pi/zhao/mc2792/lmtag/IPF_LMTAG/IPF-IND-Fingerlin_LC-All/Local/chr1.txt```

6. FDR control on significant regions using ```ipf_allregion.R```

## Step 2: Run LMTAG

1. create region file(from ipf_allregion.R). create folder to storge the output of SUPERGNOVA on selected traits 

2. create the folder that contains all phenotype-pairs results from SUPERGNOVA

2. run LMTAG pipline (/ysm-gpfs/pi/zhao/mc2792/lmtag/code)

Format: sh pipeline.sh trait region data output. example: 

```
sh pipeline.sh IPF /ysm-gpfs/pi/zhao/mc2792/lmtag/example/ipf_region_sig.txt /ysm-gpfs/pi/zhao/mc2792/lmtag/example/data /ysm-gpfs/pi/zhao/mc2792/lmtag/output2
```

3. pipeline.sh will generate the output folders and commends through main.R, then run step1.1R, step1.2R and step2.R step by step.

```/ysm-gpfs/pi/zhao/mc2792/lmtag/code/lmtag.sh``` contains the commends for execute the main LMTAG function in batch.

- lmtag.R: no filtering
- lmtag_remove.R: simply remove SNPs if W is invertible
- lmtag_ginv_noadap.R: used the same cutoff for all regions
- lmtag_ginv.R: adaptively select cutoff for different regions
- lmtag_ginv_noadap_snp.R：filtering on SNPs


## Step 3: Post GWAS analysis

### Partitioned heritability

code path: ```/ysm-gpfs/pi/zhao/yz738/Annot_MTAG/EnrichMent/code```

IMPACT annotation path: ```/ysm-gpfs/pi/zhao/yz738/software/IMPACT/IMPACT707/Annotations``` (impact707_prs_St.xlsx)

1. generate annotation:  ```gen_annot.R```

- Conditioned on baseline annotations(all annotations before genoskyline)
- IMPACT Annotation 1-707; only use EUR_; value from 0 to 1 ```IMPACT707_EUR_chr1.annot.gz```
- Create annotations for baseline annotations + each annotation of MAFK ( create one folder for each MAFK annotation)
- generated annotation: google drive and  ```/home/mc2792/scratch60/ph/annotation```

2. ldsc: ```cal_ldsc.sh```

- use the annotation generated from the previous step to calculate the ld score. make sure the output of ldsc is in the same folder of the annotation

3. partitioned heritability: ```ph.sh```

- weight and eur_chr1_snp/maf5 files: use yiliang's path




