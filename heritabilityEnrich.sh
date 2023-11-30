#!/bin/bash

### generate sumstats from OA. Tutorial command workflow is not right. The following works on 18 FEB 2022. 
### OA.txt is extracted from https://msk.hugeamp.org/downloads.html. Must include: SNP,CHR,POS,A1,A2,REF,CHOR,Beta(or OR or Z, but can only be one column otherwise error 'Too many signed sumstat columns.),P,N
### original GWAS file, only include the SNP position. Use Edirect API in NCBI to retreive relative refSNP. This procedure refer to "getSNPoa.sh".

#    python munge_sumstats.py \
#   --sumstats OAGWAS4.txt \
#    --merge-alleles w_hm3.snplist \
#    --out OA

### prepare for regression in cts mode   
filterTypeList=("FailBasic" "FailStrict")

for diseaseGroup in 0 1 100
do
for filterType in "${filterTypeList[@]}"
do 
### create a folder endoEnrich_ldscores_i for cts model heritability enrichment test
mkdir /Users/ydeng/Documents/MetaIn1.1/tempOut4In/ldsc/DiseaseGroup"$diseaseGroup"/"$filterType"/endoEnrich_ldscores

### move required files to endoEnrich_ldscores
mv /Users/ydeng/Documents/MetaIn1.1/tempOut4In/ldsc/DiseaseGroup"$diseaseGroup"/"$filterType"/*.annot.gz /Users/ydeng/Documents/MetaIn1.1/tempOut4In/ldsc/DiseaseGroup"$diseaseGroup"/"$filterType"/endoEnrich_ldscores/

mv /Users/ydeng/Documents/MetaIn1.1/tempOut4In/ldsc/DiseaseGroup$diseaseGroup/$filterType/*l2* /Users/ydeng/Documents/MetaIn1.1/tempOut4In/ldsc/DiseaseGroup$diseaseGroup/$filterType/endoEnrich_ldscores/

cp /Users/ydeng/Documents/MetaIn1.1/tempOut4In/ldsc/DiseaseGroup$diseaseGroup/$filterType/Control.GeneSet /Users/ydeng/Documents/MetaIn1.1/tempOut4In/ldsc/DiseaseGroup$diseaseGroup/$filterType/endoEnrich_ldscores/

### cts mode
### OAendo.ldcts is automatically generated from R. 
python ldsc.py \
--h2-cts OA.sumstats \
--ref-ld-chr 1000G_EUR_Phase3_baseline/baseline. \
--out /Users/ydeng/Documents/MetaOut1.3/DiseaseGroup"$diseaseGroup"/"$filterType"/DiseaseGroup"$diseaseGroup"."$filterType".EndoEnrich \
--ref-ld-chr-cts /Users/ydeng/Documents/MetaIn1.1/tempOut4In/ldsc/DiseaseGroup"$diseaseGroup"/"$filterType"/OAendo.ldcts \
--w-ld-chr weights_hm3_no_hla/weights.

echo "$filterType finished"
done

echo "Disease Group $diseaseGroup finished"
done
