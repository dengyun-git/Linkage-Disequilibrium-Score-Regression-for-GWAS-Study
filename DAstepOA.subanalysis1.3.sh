#!bin/bash

### run R script subanalysis1.3.R. This R script perform the major functionality analysis, and also generate gene annotation sets for ldsc heritability enrichment analysis.

Rscript /Users/ydeng/Documents/QCstepOA/DA_STEpUP_202111/DAcode2021.11/PrimaryAnalysis/subanalysis1.3.R

### at the same time run bash to retreive refSNP ID for all the snp locations in GWAS file 
#bash /Users/ydeng/Documents/QCstepOA/DA_STEpUP_202111/MetaIn1.1/tempOut4In/ldsc/getSNPoa.sh 
#wait


### when all the previous scripts run successfully, continuelly run the following two scripts concurrently.
### generate l2 scores for control genesets and DE genesets. 
bash l2forControl.sh &
bash l2forDE.sh
wait

### when all the l2 scores run sucessfully, run heritability enrichment test.
bash heritabilityEnrich.sh
  

