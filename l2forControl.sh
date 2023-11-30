### After preparing ensemble IDs, position on chromosomes from R. Run this l2forControl.sh to calculate l2 scores for control genesets.

### get the ldsc working environment ready. Install Anaconda. 
#conda env create --file environment.yml
#curl https://bootstrap.pypa.io/get-pip.py | python3 
#pip install pandas
#cd /Users/ydeng/Documents/QCstepOA/DA_STEpUP_202111/MetaIn1.1/tempOut4In/ldsc
#source activate ldsc
# pay attention, the arguments for files which not direct in the ldsc folder must include the full path otherwise error.

filterTypeList=("FailBasic" "FailStrict")

for diseaseGroup in 0 1 100
do
for filterType in "${filterTypeList[@]}"
do
for chrN in {1..22}
do
### generate annotation files for control proteins: Controlensemble.*.annot.gz annoys file.
python make_annot.py \
--gene-set-file  /Users/ydeng/Documents/QCstepOA/DA_STEpUP_202111/MetaIn1.1/tempOut4In/ldsc/DiseaseGroup"$diseaseGroup"/"$filterType"/Control.GeneSet \
--gene-coord-file /Users/ydeng/Documents/QCstepOA/DA_STEpUP_202111/MetaIn1.1/tempOut4In/ldsc/DiseaseGroup"$diseaseGroup"/"$filterType"/Control.Location.txt \
--windowsize 100000 \
--bimfile 1000G.EUR.QC.$chrN.bim \
--annot-file /Users/ydeng/Documents/QCstepOA/DA_STEpUP_202111/MetaIn1.1/tempOut4In/ldsc/DiseaseGroup"$diseaseGroup"/"$filterType"/Control."$chrN".annot.gz
    

### Compute partitioned LD scores with an  Controlensemble.*.annot.gz annote file. 	
python ldsc.py \
--out /Users/ydeng/Documents/QCstepOA/DA_STEpUP_202111/MetaIn1.1/tempOut4In/ldsc/DiseaseGroup"$diseaseGroup"/"$filterType"/Control."$chrN" \
--bfile 1000G.EUR.QC.$chrN \
--l2 \
--ld-wind-cm 1 \
--print-snps hm.$chrN.snp \
--annot /Users/ydeng/Documents/QCstepOA/DA_STEpUP_202111/MetaIn1.1/tempOut4In/ldsc/DiseaseGroup"$diseaseGroup"/"$filterType"/Control."$chrN".annot.gz \
--thin-annot

echo "chromosome $chrN finished"
done
echo "$filterType finished"
done
echo "Disease type $diseaseGroup finished"
done
