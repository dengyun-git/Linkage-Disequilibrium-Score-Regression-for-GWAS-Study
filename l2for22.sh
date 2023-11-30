### After preparing gene sets (Ensemble IDs) from R. Run ldsc in bash.
#!/bin/bash

### get the ldsc working environment ready. Install Anaconda. 
# conda env create --file environment.yml
# curl https://bootstrap.pypa.io/get-pip.py | python3 
# pip install pandas
#    source activate ldsc

### begin ldsc analysis.

    for chrN in {1..22}
    do
    ### generate annotation files for endo protein signatures: EDensemble.annot.gz
    python make_annot.py \
    --gene-set-file DEensembleName.GeneSet \
    --gene-coord-file DEensembleLocation.txt \
    --windowsize 100000 \
    --bimfile 1000G.EUR.QC.$chrN.bim \
    --annot-file DEensemble.$chrN.annot.gz
    
    ### Computing partitioned LD scores with an DEensemble.*.annot.gz file.
    ### strange, the command order can not work as their tutorial.The following command order works on 18 FEB 2022.
    python ldsc.py \
    --out DEendo.$chrN \
    --bfile 1000G.EUR.QC.$chrN \
    --l2 \
    --ld-wind-cm 1 \
    --print-snps hm.$chrN.snp \
    --annot DEensemble.$chrN.annot.gz \
    --thin-annot                         
    
    ### generate annotation files for control proteins: Controlensemble.*.annot.gz annoys file.
    python make_annot.py \
    --gene-set-file  Control.GeneSet \
    --gene-coord-file  ControlensembleLocation.txt \
    --windowsize 100000 \
    --bimfile 1000G.EUR.QC.$chrN.bim \
    --annot-file Controlensemble.$chrN.annot.gz
    
    ### Compute partitioned LD scores with an  Controlensemble.*.annot.gz annote file. 	
    python ldsc.py \
    --out Control.$chrN \
    --bfile 1000G.EUR.QC.$chrN \
    --l2 \
    --ld-wind-cm 1 \
    --print-snps hm.$chrN.snp \
    --annot Controlensemble.$chrN.annot.gz \
    --thin-annot
    done
    
    ### generate sumstats from OA. Tutorial command workflow is not right. The following works on 18 FEB 2022. 
    ### OA.txt is extracted from https://msk.hugeamp.org/downloads.html. Must include: SNP,CHR,POS,A1,A2,REF,CHOR,Beta(or OR or Z, but can only be one column otherwise error 'Too many signed sumstat columns.),P,N
    ### original GWAS file, only include the SNP position. Use Edirect API in NCBI to retreive relative refSNP. This procedure refer to "getSNPoa.sh".
    python munge_sumstats.py \
    --sumstats OAGWAS3.txt \
    --merge-alleles w_hm3.snplist \
    --out OA
    
    ### create a folder endoEnrich_ldscores for cts model heritability enrichment test
    mkdir endoEnrich_ldscores
    ### move required files to endoEnrich
    mv DEensemble.*.annot.gz DEendo.*.l2.ldscore.gz DEendo.*.l2.M DEendo.*.l2.M_5_50 Controlensemble.*.annot.gz Control.*.l2.ldscore.gz Control.*.l2.M Control.*.l2.M_5_50 endoEnrich_ldscores/
    cp Control.GeneSet endoEnrich_ldscores/

    ### cts mode
    ### OAendo.ldcts should be pre-edited as: EndoType	endoEnrich_ldscores/DEendo.,endoEnrich_ldscores/Control. 
    python ldsc.py \
    --h2-cts OA.sumstats.gz \
    --ref-ld-chr 1000G_EUR_Phase3_baseline/baseline. \
    --out endo1 \
    --ref-ld-chr-cts OAendo.ldcts \
    --w-ld-chr weights_hm3_no_hla/weights.
   
