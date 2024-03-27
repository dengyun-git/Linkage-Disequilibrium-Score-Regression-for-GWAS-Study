# Subanalysis 1.3 Bioinformatic characterization for proteins signatures distinguishing clusters
# Date: 16 MAR 2022
# Version: 0.7
# Author(s): name = Dr. Yun Deng
# 
# Description: The purpose of this R script is to analyse bioinformatic characteristics of protein signatures for each entodype 
# inputs & outputs file lists refer to log1.3

### Load relevant libraries
library(ggplot2)
library(factoextra)
library(magrittr)
library(WGCNA)
library(SomaDataIO)
library(readxl) 
library(fgsea)
library(org.Hs.eg.db)
library(scuttle)
# library(EWCE)
library(haplo.stats)
library(QuaternaryProd)
library(enrichplot)
library(cowplot)
library(ggridges)
library(ggnewscale)
library(biomaRt)
library(grid)
library(gridExtra)
library(gplots)
library(data.table)
options(stringsAsFactors = FALSE)
library(colorspace)

### set file path
myPathIn1  <- "/Users/ydeng/Documents//DA__202111/MetaIn1.1/_DAG_rel001_draft1/"  ### location of downloaded meta data from OneDrive
myPathIn2 <- "/Users/ydeng/Documents//DA__202111/MetaIn1.1/DatabaseResource/"         ### location of downloaded enrichment resources from online databases
myPathIn3  <- "/Users/ydeng/Documents//DA__202111/MetaOut1.2/normalised_combat/"  
myPathIn4  <- "/Users/ydeng/Documents//DA__202111/MetaOut1.1/normalised_combat/"       ### location of files required as input in this script
myPathIn5 <- "/Users/ydeng/Documents//DA__202111/MetaIn1.1/tempOut4In/"   
### location of files required as input in this script
myCodeIn  <- "/Users/ydeng/Documents//DA__202111/DAcode2021.11/"                       ### location of downloaded R scripts from Github

source(paste0(myCodeIn,"PrimaryAnalysis/call1.1.R")) 
source(paste0(myCodeIn,"PrimaryAnalysis/call1.3.R")) 

myPathOut <- "/Users/ydeng/Documents//MetaOut1.3/"

### create a log file to record the "code version + running date + input files + output files"
log1.3  = file(paste0(myPathOut,"1.3.log"),"append")
### write to the log file to trace the code running and developing history 
writeLines(paste0("Running time ",Sys.Date(),".\n"),log1.3)
writeLines(paste0("Code version: 0.43 \n"),log1.3)

DiseaseGroupList = c(0,1,100)                   ### disease groups: 1 = OA | 0 = Injury | 2 = Control | 3 = Inflammatory arthritis control | 100 = OA + Injury
filterTypeList = c("FailBasic","FailStrict")    ### protein filter in QC procedure
RFUstart = "CRYBB2.10000.28"                    ## the first protein in the Somologic data file


### read in previous files, download required files and process them
###1. read in raw protein meta data, containing information from somalogic adat file.##################
ProColMeta <- read.csv(paste0(myPathIn1,"somascan/csv/ProColMeta.txt"),sep=";")
writeLines(paste0("Input file:",paste0(myPathIn1,"somascan/csv/ProColMeta.txt"),"\n"),log1.3)

### remove observations with empty "EntrezGeneID"
# ProColMeta <- ProColMeta[-which(ProColMeta$EntrezGeneID==""),] 
### For the purpose of extract information for multiple genes encoding complex proteins, tidy up the EntrezGeneID, EntrezGeneSymbol as well as UniProt which have multiple spaces
for(Counter in 1:nrow(ProColMeta)){
  ProColMeta$UniProt[Counter] <- sub("\\s+","",ProColMeta$UniProt[Counter]) ### remove space only leave "," as delimiter
  ProColMeta$EntrezGeneID[Counter] <- sub("\\s+"," ",ProColMeta$EntrezGeneID[Counter]) ### merge multiple spaces as one single space
  ProColMeta$EntrezGeneSymbol[Counter] <- sub("\\s+"," ",ProColMeta$EntrezGeneSymbol[Counter])
}

### find proteins with multiple EntrezGenes
dupRowIdA <- grep(" ",ProColMeta$EntrezGeneID)
### a data frame showing proteins with multiple EntrezGenes. proComplexRef is used to replace element in the reference gene sets. 
proComplexRef = ProColMeta[dupRowIdA,c("SomaId","EntrezGeneID","EntrezGeneSymbol","UniProt")]


### read in protein expression profile (code reference subanalysis 1.1)
### Read in proteomic profile according to user defined normalization methods
SomaFrame <- read.csv(paste0(myPathIn1,"somascan/csv/SomaScan_rel1_normalised_combat.csv")) ### read in proteomics expression data which underwent normalisation + combatted correction
rownames(SomaFrame) <- SomaFrame$X
writeLines(paste0("Input file:",paste0(myPathIn1,"somascan/csv/SomaScan_rel1_normalised_combat.csv"),"\n"),log1.3)

### Read in clinical and QC data (& give the dataframes their row names (STEp-UP OA IDs))
QCFrame <- read.csv(paste0(myPathIn1,'clinical/discovery_QApheno_1.csv'))
rownames(QCFrame) = QCFrame$sf_iknee_sample_id_number
ClinicFrame1 <- read.csv(paste0(myPathIn1,'clinical/discovery_DAPpheno_1.csv'))
rownames(ClinicFrame1) =ClinicFrame1$sf_iknee_sample_id_number
ClinicFrame2 <- read.csv(paste0(myPathIn1,'clinical/discovery_DAPpheno_2.csv'))
rownames(ClinicFrame2) =ClinicFrame2$sf_iknee_sample_id_number
### Combine clinical and QC dataframes and to ensure the row order matches the Somalogic dataframe.
CombinedFrame<- getMyCombinedF(QCFrame,ClinicFrame1,ClinicFrame2,SomaFrame)
### adjust proteomic expression scale because combat correction has applied log transformation
CombinedFrame[,which(colnames(CombinedFrame)==RFUstart):ncol(CombinedFrame)] <- CombinedFrame[,which(colnames(CombinedFrame)==RFUstart):ncol(CombinedFrame)] %>% exp()
writeLines(paste0("Input file:",paste0(myPathIn1,'clinical/discovery_QApheno_1.csv'),"\n"),log1.3)

### Load filters 
sampFilter <- read.csv(paste0(myPathIn1,'somascan/QC/sample_filters.txt'),sep='\t')
protFilter <- read.csv(paste0(myPathIn1,'somascan/QC/protein_filters.txt'),sep='\t')
writeLines(paste0("Input files:",paste0(myPathIn1,'somascan/QC/sample_filters.txt'),"\n",
                  paste0(myPathIn1,'somascan/QC/protein_filters.txt'),"\n"),log1.3)



### 2. downloaded csv from https://www.proteinatlas.org/about/download ########################
SubLocation <- read.csv(paste0(myPathIn2,"subcellular_location.tsv"),sep="\t") ### consider whether to take into account additional location column

LocationList <- diffProCombineOnline(SubLocation,"Main.location")
SubLocationDat = LocationList[[1]]
SubLocationType = LocationList[[2]]
SublocationSet = LocationList[[3]]

### supplement exosome database into SublocationSet list 
Exosomoe = read.table(paste0(myPathIn2,"exosome.tsv"),header=TRUE,sep="\t")$GeneSymbol
### add exosome to the sublocation gene sets
SublocationSet[["Exosomoe"]] = Exosomoe
### protein complexies
SublocationSet2 <- replaceGenes(SublocationSet,proComplexRef,"EntrezGeneSymbol") 


### broader cell location type enrichment test
Cytoplasm <- scan(paste0(myPathIn2,"Cytoplasm.txt"),what ="character",quiet=TRUE)
Nucleus <- scan(paste0(myPathIn2,"Nucleus.txt"),what ="character",quiet=TRUE)
Endomembrane <- scan(paste0(myPathIn2,"Endomembrane.txt"),what ="character",quiet=TRUE)
secreted <- read.table(paste0(myPathIn2,"secreted_protein.tsv"),sep="\t",colClasses = c("NULL","character",rep("NULL", 10)), header = TRUE)

### for each gene map their detailed cell location to broader cell location types. 
for (indexCounter in 1:nrow(SubLocationDat)){
  x = SubLocationDat$Main.location[indexCounter]
  mainLoc = strsplit(x,";")[[1]]
  if(any(Cytoplasm %in% mainLoc)){newLocation = "Cytoplasm"
  }else if(any(Nucleus %in% mainLoc)){newLocation = "Nucleus"
  }else if(any(Endomembrane %in% mainLoc)){newLocation = "Endomembrane"
  }else{newLocation = "secreted"}
  SubLocationDat$Broad.location[indexCounter]= newLocation
}

BroadLocationList <- diffProCombineOnline(SubLocationDat,"Broad.location")
SubLocationBroad = BroadLocationList[[1]]
SubLocationBroadType = BroadLocationList[[2]]
SubLocationBroadSet = BroadLocationList[[3]]
### protein complexies
SubLocationBroadSet2 <- replaceGenes(SubLocationBroadSet,proComplexRef,"EntrezGeneSymbol")


### 3. retreive cell tissue tyeps information ########################
cell.tissue <- read.csv(paste0(myPathIn2,"normal_tissue.tsv"),sep="\t")
cell.tissue[,"Tissue"] <-  gsub('[[:digit:]]', '',cell.tissue[,"Tissue"]) ### clear unused suffix numbers in "Tissue" column
cell.tissue <- cell.tissue[-which(cell.tissue$Tissue == "N/A" | cell.tissue$Cell.type=="N/A"),] ### clear N/A rows

### matching sigPrTableUni to tissueType
TissueList <- diffProCombineOnline(cell.tissue,"Tissue")
TissueDat = TissueList[[1]]
TissueType = TissueList[[2]]
TissueSet = TissueList[[3]]
### protein complexies
TissueSet2 <- replaceGenes(TissueSet,proComplexRef,"EntrezGeneSymbol")

CellList <- diffProCombineOnline(cell.tissue,"Cell.type")
CellDat = CellList[[1]]
CellType = CellList[[2]]
CellSet = CellList[[3]]
CellSet2 <- replaceGenes(CellSet,proComplexRef,"EntrezGeneSymbol")


### 4. Apopsis Necorosis analysis ########################
### use biomaRt package to find homologous genes between mouse and human
### Map mouse EntrezGeneID from GI using uniProt Retrieve/ID mapping tool. Refere to "ConvertGI.sh"
### 21FEB 2022, also try a more standard tool "NCBI Entrez Programming Utilities" to convert GI to Entrez. https://www.ncbi.nlm.nih.gov/books/NBK25498/#chapter3.Application_1_Converting_GI_num
### Besides previous methods, also can GI -> ENTREZ ID, use EDirect NCBI API.

### find which Mart and which dataset to use:
# biomaRt::listMarts()
# myMart = biomaRt::useMart("ENSEMBL_MART_ENSEMBL")
# biomaRt::listDatasets(myMart)
# biomaRt::searchDatasets(mart = myMart, pattern = "hsapiens")
### find which attributes and filters to use for function getLDS (considering our case is: EntrezGeneID is the most convenient term to use):
### listAttributes(humGenes)[grep("entre",biomaRt::listAttributes(humGenes)[,1]),] 
### listFilters(humGenes)[grep("entre",biomaRt::listFilters(humGenes)[,1]),] 

humGenes = biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
musGenes = biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")

### Generate Necrosis gene sets, Apoptosis gene sets. "ApoptosisEntrezIDmouse.csv" and "NecrosisEntrezIDmouse.csv" contains human homoglous genes mapped from mouse genes.
NecGeneSet <- getMyHumGene(musGenes,humGenes,paste0(myPathIn2,"/NecrosisApoptosisAlarming/NecrosisUniProID.csv"))
ApopGeneSet <- getMyHumGene(musGenes,humGenes,paste0(myPathIn2,"/NecrosisApoptosisAlarming/ApoptosisUniProID.csv"))
### Generate Alarming gene sets.
AlarGeneName <- unique(scan(paste0(myPathIn2,"/NecrosisApoptosisAlarming/Alarmings.csv"),sep=",", what ="character",quiet=TRUE))
AlarGeneSet <- getBM(attributes = c('entrezgene_id'),filters = 'entrezgene_accession', values = AlarGeneName , mart = humGenes)
### use myHumGene to generate the Necrosis|Apoptosis|Alarming gene sets.
NecApopAlarGeneSetPre = list("Necrosis"=as.character(NecGeneSet),"Apoptosis"=as.character(ApopGeneSet),"Alarming"=unlist(AlarGeneSet))
### tackle with protein complexies
NecApopAlarGeneSet <- replaceGenes(NecApopAlarGeneSetPre,proComplexRef,"EntrezGeneID")


### 5. read in raw scRNA expression matrix, using QC matrix provided by the paper to filter poor quaility cells and genes.Obtain scRNAexp as the count matrix for downstream analysis.########################
seqMetaF <- read.table(gzfile(paste0(myPathIn2,"celseq_meta.tsv.725591.gz")),sep='\t',header=TRUE)
remainCell <- seqMetaF$cell_name[which(seqMetaF$percent_mt_molecules <= 0.25 & seqMetaF$genes_detected > 1000 &seqMetaF$disease=="OA")] #discarded cells with fewer than 1,000 genes detected with at least one fragment & cells that had more than 25% of molecules coming from mitochondrial genes & cells for OA patients. 
scRNAexpF <- read.table(gzfile(paste0(myPathIn2,"celseq_matrix_ru10_molecules.tsv.725585.gz")),sep="\t",header=TRUE)
remainGene <- apply(scRNAexpF[,-1],1,function(x){length(which(!is.na(x)))>10}) #discarded genes that had nonzero expression in fewer than 10 cells.
scRNAexp <- scRNAexpF[remainGene,remainCell] ### up to here, simple QC done. 
rownames(scRNAexp) <- scRNAexpF[,"gene"][remainGene]
scRNAexp[is.na(scRNAexp)]=0

### map cells to cell types
cells = colnames(scRNAexp)
cellType <- seqMetaF$type[sapply(cells,function(x){which(seqMetaF$cell_name==x)})]
colnames(scRNAexp) <- cellType
cellType2 = unique(cellType)
###new scRNAexp based on cellType
scRNAexp2 = matrix(NA, ncol=length(cellType2),nrow=nrow(scRNAexp))
rownames(scRNAexp2) = rownames(scRNAexp)
for(scCounter in 1:length(cellType2)){
  scRNAexp2[,scCounter] = rowSums(scRNAexp[,which(colnames(scRNAexp)==cellType2[scCounter])])
}
colnames(scRNAexp2) = cellType2

OA_rna = list(exp=as.matrix(scRNAexp2),annot=list(level1class=cellType2))

### Generate cell type data 
annotLevels <- list(level1class=OA_rna$annot$level1class)
fNamesOArna <- EWCE::generate_celltype_data(exp = OA_rna$exp,annotLevels = annotLevels,groupName = "OA") 
file.copy(fNamesOArna, myPathIn2)
load(paste0(myPathIn2,"CellTypeData_OA.rda"))
file.remove(paste0(myPathIn2,"CellTypeData_OA.rda"))  ### clear the file for the next round copy, otherwise file.copy will fail

### get the specificity frame
specificityFrame <- ctd[[1]]$specificity
genes <- rownames(specificityFrame)

### generate cell type specific marker gene sets using top 10% specificity genes for each cell type
scRNAgeneSets = list()
for(cellCounter in 1:length(cellType2)){
  geneInd <- order(-specificityFrame[,cellCounter])[1:floor(0.1*length(genes))]
  scRNAgeneSets[[cellCounter]] =  genes[geneInd]
  names(scRNAgeneSets)[cellCounter] = cellType2[cellCounter]
}

### exchange cell names to cell types, each cell type include gene names
scRNAgeneCellSets = tapply(scRNAgeneSets,names(scRNAgeneSets),function(x) unique(unlist(x)))
### tackle with protein complexies
scRNAgeneCellSets2 = replaceGenes(scRNAgeneCellSets,proComplexRef,"EntrezGeneSymbol")

### 6. GWAS file preparation. For the convenience, use both bash and R.
### bash: combine the chromosom position information from bim file
# cat 1000G.EUR.QC.bim | awk '{print $1":"$4"_"$5"_"$6"\t"$2}' > snpMap.txt
# ### bash:extract the chromosome position information from gwas file
# cat KP.Format.GO.FILTER.GW.AllOA.FULL.09052019.txt | cut -f 1 | sed -e '1d' > GWASposition.txt
# ### R: match the snp based on position in gwas file
# snpMap <- fread("/Users/ydeng/Documents//DA__202111/MetaIn1.1/tempOut4In/ldsc/snpMap.txt",head=FALSE)
# GWASposition <- fread("/Users/ydeng/Documents//DA__202111/MetaIn1.1/tempOut4In/ldsc/GWASposition.txt",head=FALSE)
# SNP <- unlist(sapply(GWASposition,function(x){snpMap[which(snpMap[,1] %in% x),2]}))
# colnames(SNP)="SNP"
# fwrite(SNP,"/Users/ydeng/Documents//DA__202111/MetaIn1.1/tempOut4In/ldsc/SNP.txt",sep="\n",col.names=TRUE)
# ### bash: combine the mapped refSNP id to "KP.Format.GO.FILTER.GW.AllOA.FULL.09052019.txt"
# cut -f2- KP.Format.GO.FILTER.GW.AllOA.FULL.09052019.txt > temptext.txt 
# paste SNP.txt temptext.txt > OAGWAS.txt


writeLines(paste0("Input files:",paste0(myPathIn2,"subcellular_location.tsv"),"\n",
                  paste0(myPathIn2,"exosome.tsv"),"\n",
                  paste0(myPathIn2,"Cytoplasm.txt"),"\n",
                  paste0(myPathIn2,"Nucleus.txt"),"\n",
                  paste0(myPathIn2,"Endomembrane.txt"),"\n",
                  paste0(myPathIn2,"secreted_protein.tsv"),"\n",
                  paste0(myPathIn2,"normal_tissue.tsv"),"\n",
                  paste0(myPathIn2,"/NecrosisApoptosisAlarming/NecrosisUniProID.csv"),"\n",
                  paste0(myPathIn2,"/NecrosisApoptosisAlarming/ApoptosisUniProID.csv"), "\n",
                  paste0(myPathIn2,"/NecrosisApoptosisAlarming/Alarmings.csv"),"\n",
                  paste0(myPathIn2,"celseq_meta.tsv.725591.gz"),"\n",
                  paste0(myPathIn2,"h.all.v7.5.1.entrez.gmt.txt"),"\n",
                  paste0(myPathIn2,"c2.cp.kegg.v7.5.1.entrez.gmt.txt"),"\n",
                  paste0(myPathIn2,"c5.go.v7.5.1.entrez.gmt.txt")),log1.3)


### for loop begin#############################
for(diseaseGroup in DiseaseGroupList){
  
  ### create a folder per disease group within corresponding normalization folders
  dir.create(paste0(myPathOut,"DiseaseGroup",diseaseGroup))
  ### set the gwas ldsc directory to take control gene sets file
  dir.create(paste0(myPathIn5,"ldsc/DiseaseGroup",diseaseGroup))
  
  for(filterType in filterTypeList){
    
    ### create folder per filter type within corresponding "normalisation/diseaseGroup/"
    dir.create(paste0(myPathOut,"DiseaseGroup",diseaseGroup,"/",filterType))
    outF1 = paste0(myPathOut,"DiseaseGroup",diseaseGroup,"/",filterType,"/")
    
    ### generate sigPrTable: contain differentially expression statistics and protein meta table information such as keys in different database 
    ProTable <- read.csv(paste0(myPathIn3,"DiseaseGroup",diseaseGroup,"/",filterType,"/AbundanceTest.csv"),sep=",")
    writeLines(paste0("Input file:",paste0(myPathIn3,"DiseaseGroup",diseaseGroup,"/",filterType,"/AbundanceTest.csv"),"\n"),log1.3)
    
    ### combine differentially expression test statistics with corresponding protein meta data 
    sigProColMeta <- ProColMeta[rownames(ProTable),]  ### ProColMeta does not has the same row numbers as ProTable, because ProColMeta include unfiltered proteins
    sigPrTableR <- cbind(ProTable,sigProColMeta)
    ### remove NA rows
    sigPrTableR <- sigPrTableR[-which(is.na(sigPrTableR$EntrezGeneID)),]
    
    ### enrichment analysis "duplicate gene names, fgsea may produce unexpected results", we keep the set of unique entries in "sigPrTableUni".
    sigPrTableUni = sigPrTableR[sapply(unique(sigPrTableR[,"EntrezGeneID"]),function(x){which(sigPrTableR[,"EntrezGeneID"]==x)[1]}),]
    
    ### we also need protein expression profile in the downstream analysis
    ### apply filters: sample filters and protein filters
    if(filterType=="FailBasic"){
      ### for basic fitlers, filter samples
      SampExclude <- sampFilter$Sample[sampFilter$inj_filter_status==filterType|sampFilter$oa_filter_status==filterType]
      
      ### for basic filters, filter proteins
      ProtExclude <- protFilter$Protein[protFilter$inj_filter_status==filterType|protFilter$oa_filter_status==filterType]
    }else{
      ### for strict filters, filter samples
      SampExclude <- sampFilter$Sample[c(grep("Fail",sampFilter$inj_filter_status),grep("Fail",sampFilter$oa_filter_status))]
      
      ### for strict filters, filter proteins
      ProtExclude <- protFilter$Protein[c(grep("Fail",protFilter$inj_filter_status),grep("Fail",protFilter$oa_filter_status))]
    }
    
    ### obtain RUF frame after filtering: FrameFil
    FrameFil <- CombinedFrame[!(CombinedFrame$sf_iknee_sample_id_number%in%SampExclude),!(colnames(CombinedFrame)%in%ProtExclude)]
    
    ### filtered proteomic expression frame further down to baseline sample, select current disease group
    if(diseaseGroup==100){ProExpReady <- FrameFil[FrameFil$baseline==1 & (FrameFil$sf_iknee_qc_group == 0|FrameFil$sf_iknee_qc_group == 1),]
    }else{ProExpReady <- FrameFil[FrameFil$baseline==1 & FrameFil$sf_iknee_qc_group == diseaseGroup,]}
    
    ### also read in endotypes
    ### read in endotypes
    heatEndoF <- read.csv(file=paste0(myPathIn4,"DiseaseGroup",diseaseGroup,"/",filterType,"/normalised_combatDis",diseaseGroup,filterType,"PCAumapEndo.csv"))
    ### matching heatEndoF with ProExpReady
    EndoAll <- heatEndoF[rownames(ProExpReady),which(colnames(heatEndoF)=="EndoLabel")] ### Make sure the rownames are matching between heatEndo and topkExpFrame
    writeLines(paste0("input file:",paste0(myPathIn4,"DiseaseGroup",diseaseGroup,"/",filterType,"/normalised_combatDis",diseaseGroup,filterType,"PCAumapEndo.csv"),"\n"),log1.3)
    
    ### calculate fold change between endotypes for each protein
    fcFrame <- as.matrix(colMeans(ProExpReady[which(EndoAll==1),which(colnames(ProExpReady)==RFUstart):ncol(ProExpReady)])/colMeans(ProExpReady[which(EndoAll!=1),which(colnames(ProExpReady)==RFUstart):ncol(ProExpReady)]),ncol=1)
    ### combine fcFrame to sigPrTableUni
    sigPrTableUni <- cbind(sigPrTableUni,fcFrame[rownames(sigPrTableUni),])
    colnames(sigPrTableUni)[ncol(sigPrTableUni)]="fc"
    fc = sigPrTableUni$fc
   
    ### outside pCol loop, for GWAS, convert entrezGene name to ENSEMBLE GENE ID. differentially expressed "DEensembleName.txt", control genesets in "ControlensembleName.txt"
    ### control set
    HUGOcontrol <- unique(getBM(attributes=c("hgnc_symbol"),filters = "entrezgene_accession", values=sigPrTableUni$EntrezGeneSymbol, mart= humGenes)) 
    ControlEnsemble <- getBM(attributes=c("ensembl_gene_id","chromosome_name","start_position", "end_position"),filters = "hgnc_symbol", values=HUGOcontrol, mart= humGenes)
    ControlEnsemble2 <- ControlEnsemble[which(!grepl("_|X",ControlEnsemble$chromosome_name)),]
    colnames(ControlEnsemble2) = c("GENE","CHR","START","END")
    ControlEnsemble2$CHR <- paste0("ch",ControlEnsemble2$CHR)
    ControlEnsembleName <- as.matrix(ControlEnsemble2[,1],ncol=1)
    
    dir.create(paste0(myPathIn5,"ldsc/DiseaseGroup",diseaseGroup,"/",filterType))
    outF2 = paste0(myPathIn5,"ldsc/DiseaseGroup",diseaseGroup,"/",filterType)
    write.table(ControlEnsembleName, paste0(outF2,"/Control.GeneSet"), row.names = FALSE, col.names=FALSE, quote=FALSE)
    write.table(ControlEnsemble2,paste0(outF2,"/Control.Location.txt"),row.names = FALSE,col.names=TRUE,quote=FALSE)
    writeLines(paste0("Onput file:",paste0(outF2,"/Control.GeneSet"),"\n"),log1.3)
    writeLines(paste0("Onput file:",paste0(outF2,"/Control.Location.txt"),"\n"),log1.3)
    
    ### run bash script to generate LD scores for control gene sets.
    
    ### required data frame get ready for enrichment analysis########################
    ### gene list will be ranked according to different pvalues derived from different regression methods. Namely, columns in sigPrTableUni
    pEnd = grep("Pvalue",colnames(sigPrTableUni))
    for(pColum in pEnd){
      
      rankPro = sigPrTableUni[,pColum]
     
      if(any(is.na(rankPro))){next      ### if NA,exist such p value column
      }else{
        ### create the plot output file
        pdf(file=paste0(outF1,"Dis",diseaseGroup,filterType,colnames(sigPrTableUni)[pColum],".pdf"))
        writeLines(paste0(outF1,"Dis",diseaseGroup,filterType,colnames(sigPrTableUni)[pColum],".pdf","\n"),log1.3)
        
        
        ### add Enrichment test Gene symbols to sigPrTableUni after considering protein complexity
        ### "EnrichEntrezIDLabel" and "EnrichEntrezSymbolLabel" columns in sigPrTableUni are used as names for ranked p vlaues for enrichment test
        dupRowId <- grep(" ",sigPrTableUni$EntrezGeneID)
        names(rankPro)[dupRowId] <- sigPrTableUni$SomaId[dupRowId]
        names(rankPro)[-dupRowId] <-  sigPrTableUni$EntrezGeneID[-dupRowId]
        sigPrTableUni[,"EnrichEntrezIDLabel"] = names(rankPro)
        
        names(rankPro)[dupRowId] <- sigPrTableUni$SomaId[dupRowId]
        names(rankPro)[-dupRowId] <-  sigPrTableUni$EntrezGeneSymbol[-dupRowId]
        sigPrTableUni[,"EnrichEntrezSymbolLabel"] = names(rankPro)
        
        ### pathway analysis using MSigDB database. Names of the ranked p vluaes using combination of EntrezGeneID and SomaId(genes encoding complex proteins)##########

        ### for 2 endotypes, select the down regulated genes as endotype marker proteins.
        IDremain <- which(rankPro!=0 & fc!=0)
        if(length(unique(EndoAll))==2 & colnames(sigPrTableUni)[pColum]=="Logistic.Pvalue.Endo.2"){fc = 1/fc}
        rankPro1 <- -log(rankPro[IDremain])*sign(log(fc[IDremain]))
        names(rankPro1) <- sigPrTableUni$EnrichEntrezIDLabel[IDremain]
        ranks1 <- sort(rankPro1,decreasing=TRUE)
        
        ### pathway enrichment adopt four databases, locally downloaded as h.all.v7.5.1.entrez.gmt.txt; c2.cp.kegg.v7.5.1.entrez.gmt.txt; c5.go.v7.5.1.entrez.gmt.txt.
        padjThresh = 0.5
        ###Hall database
        HallEnrichList<- getEnrichPath(myPathIn2,"h.all.v7.5.1.entrez.gmt.txt",proComplexRef,ranks1,1.1)
        HallPath <- HallEnrichList[[1]]
        HallGeneSets <- HallEnrichList[[2]]
        HallMain <- HallEnrichList[[3]]
        topPathways2 <- HallMain[which(HallMain$padj<padjThresh), pathway]
        # pretty plot when combining these two
        # gseaplot2(HallPath, geneSetID = 1:2, pvalue_table = TRUE,color = c("#E495A5", "#86B875", "#7DB0DD"), ES_geom = "dot")
        # plotMyGseaTable("Pathway",HallGeneSets[topPathways2], ranks1, HallPath,gseaParam=1)
        fwrite(HallPath, file=paste0(outF1,"Dis",diseaseGroup,filterType,colnames(sigPrTableUni)[pColum],"HallPath.txt"), sep="\t", sep2=c("", " ", ""))
        fwrite(HallMain, file=paste0(outF1,"Dis",diseaseGroup,filterType,colnames(sigPrTableUni)[pColum],"HallMain.txt"), sep="\t", sep2=c("", " ", ""))
        writeLines(paste0("Onput file:",paste0(outF1,"Dis",diseaseGroup,filterType,colnames(sigPrTableUni)[pColum],"HallPath.txt"),"\n"),log1.3)
        
        ### KEGG database
        KEGGEnrichList<- getEnrichPath(myPathIn2,"c2.cp.kegg.v7.5.1.entrez.gmt.txt",proComplexRef,ranks1,1.1)
        KEGGPath <- KEGGEnrichList[[1]]
        KEGGGeneSets <- KEGGEnrichList[[2]]
        KEGGMain <- KEGGEnrichList[[3]]
        topPathways3 <- KEGGMain[which(KEGGMain$padj<padjThresh),pathway]
        # plotMyGseaTable("Pathway",KEGGGeneSets[topPathways3], ranks1, KEGGPath,gseaParam=0.2)
        fwrite(KEGGPath, file=paste0(outF1,"Dis",diseaseGroup,filterType,colnames(sigPrTableUni)[pColum],"KEGGPath.txt"), sep="\t", sep2=c("", " ", ""))
        fwrite(KEGGMain, file=paste0(outF1,"Dis",diseaseGroup,filterType,colnames(sigPrTableUni)[pColum],"KEGGMain.txt"), sep="\t", sep2=c("", " ", ""))
        writeLines(paste0("Onput file:",paste0(outF1,"Dis",diseaseGroup,filterType,colnames(sigPrTableUni)[pColum],"KEGGPath.txt"),"\n"),log1.3)
        
        ### GO database
        GOEnrichList<- getEnrichPath(myPathIn2,"c5.go.v7.5.1.entrez.gmt.txt",proComplexRef,ranks1,1.1)
        GOPath <- GOEnrichList[[1]]
        GOGeneSets <- GOEnrichList[[2]]
        GOMain <- GOEnrichList[[3]]
        topPathways4 <- GOMain[which(GOMain$padj<padjThresh), pathway]
        # plotMyGseaTable("Pathway",GOGeneSets[topPathways4], ranks1, GOPath,gseaParam=0.2)
        fwrite(GOPath, file=paste0(outF1,"Dis",diseaseGroup,filterType,colnames(sigPrTableUni)[pColum],"GOPath.txt"), sep="\t", sep2=c("", " ", ""))
        fwrite(GOMain, file=paste0(outF1,"Dis",diseaseGroup,filterType,colnames(sigPrTableUni)[pColum],"GOMain.txt"), sep="\t", sep2=c("", " ", ""))
        writeLines(paste0("Onput file:",paste0(outF1,"Dis",diseaseGroup,filterType,colnames(sigPrTableUni)[pColum],"GOPath.txt"),"\n"),log1.3)
        
        ### enriched pathway visualisation############
        
        ### gene map: gene membership to pathways & display pathway crosstalk. 
        ### using clusterProfiler for visualisation
        try({p22 <- CrossTalkMap(myPathIn2,"h.all.v7.5.1.entrez.gmt.txt",topPathways2,ranks1,fc)
        print(p22)})
        try({p23 <- CrossTalkMap(myPathIn2,"c2.cp.kegg.v7.5.1.entrez.gmt.txt",topPathways3,ranks1,fc)
        print(p23)})
        try({p24 <- CrossTalkMap(myPathIn2,"c5.go.v7.5.1.entrez.gmt.txt",topPathways4,ranks1,fc)
        print(p24)})
        
        
        ### heatmap for highly differentially expressed genes and their corresponding enriched pathways
        ### extract topk highly expressed gene names from sigPrTableUni, use the gene names in somalogic adat file.
        topk = 100
        topkGene <- rownames(sigPrTableUni[order(rankPro)[1:topk],]) ### when sigPrTableUni was generated it only contains the filtered proteins, so topK genes also are within ProExpReady

        ### generate dataframe showing expression level of highly expressed proteins. map topk to names in somalogic adat file,  but protein names use Entreze Gene ID (for the convenience of mapping to pathways)
        topkExpFrame <- ProExpReady[,topkGene]
        ### top gene name use entrez id with protein complixies replacements
        ### topkGene2 is used to map with gene sets; colnames(topkExpFrame) is used to display in plot
        colnames(topkExpFrame) <- sigPrTableUni[order(rankPro)[1:topk],"EnrichEntrezSymbolLabel"]
        topkGene2 = sigPrTableUni[order(rankPro)[1:topk],"EnrichEntrezIDLabel"] ### EntrezGeneID matching colnames(topkExpFrame)

        ### map top differentailly expressed genes to their enriched pathway categories
        PathColor2 <- getColBarColor(HallPath,topPathways2,HallGeneSets,topkGene2,"cyan3")
        PathColor3 <- getColBarColor(KEGGPath,topPathways3,KEGGGeneSets,topkGene2,"darkseagreen2")
        PathColor4 <- getColBarColor(GOPath,topPathways4,GOGeneSets,topkGene2,"deeppink4")
       
        ### generate clab paramerter for heatmap3
        clab=t(rbind(PathColor2,PathColor3,PathColor4))
        pathLegend=c("NA",rownames(PathColor2),rownames(PathColor3),rownames(PathColor4))
        legendColor=c("grey","cyan3","darkseagreen2","deeppink4")

        ### generate rlab parameter (representing endotypes) for heatmap3
        heatEndo <- heatEndoF[rownames(topkExpFrame),"EndoLabel"] ### Make sure the rownames are matching between heatEndo and topkExpFrame
        endoColor = colorspace::sequential_hcl(length(unique(heatEndo)), palette = "Hawaii")
        endoType = unique(heatEndo)

        ### for better visualisation, put the rows of the same endos together
        newOrder = unlist(sapply(endoType,function(x){which(heatEndo==x)}))
        heatEndo <- heatEndo[newOrder]
        row_annotation = as.matrix(t(unlist(sapply(heatEndo,function(x){endoColor[which(endoType==x)]}))))
        rownames(row_annotation)=c("entotype")

        ### adjust the expression matrix accordingly
        topkExpFrame2 = topkExpFrame[newOrder,]

        ### clustered heatmap on highly differentailly expressed topGene2, column bar color showing membership for clusters(if any enriched pathways have been found)
        heatmap.3(sweep(sweep(log(topkExpFrame2),2,colMeans(as.matrix(log(topkExpFrame2))),"-"),2,colSds(as.matrix(log(topkExpFrame2))),"/"), scale="none",
                  dendrogram="column", margins=c(5,9),Rowv=FALSE, Colv=TRUE,
                  ColSideColors=clab, RowSideColors=row_annotation, symbreaks=FALSE, ColSideColorsSize=3, RowSideColorsSize=0.5,
                  key=TRUE, symkey=FALSE, density.info="none", trace="none", KeyValueName="Normalized log Expression Level",
                  lhei=c(1,4), lwid=c(1,2.5),
                  main=paste0("Enriched pathway for top differentially expressed proteins \n based on Disease Group", diseaseGroup," ",filterType," ",colnames(sigPrTableUni)[pColum]),
                  labCol=colnames(topkExpFrame2),labRow=rownames(topkExpFrame2), cexRow=0.7,cexCol = 0.7,col=rev(heat.colors(75)))
        legend(y=0.5, x=-.07, xpd=TRUE,legend=pathLegend,fill=legendColor,
               border=FALSE, bty="n",cex=0.6,lwd=0.1,text.width=0.6) ### xpd self define location to plot legend
        # 
        # 
        ### cell sublocation enrichment analysis##########
        ### prepare ranked gene list for sublocation enrichment test
        rankPro3 <- -log(rankPro[IDremain])*sign(log(fc[IDremain]))
        names(rankPro3) <- sigPrTableUni$EnrichEntrezSymbolLabel[IDremain]
        ranks3 <- sort(rankPro3)
        
        ### subcellular location enrichment test
        fgseaLocal <- fgsea(pathways = SublocationSet2, stats=ranks3, scoreType = "std",eps = 0.0,minSize=10, maxSize=1000)
        sigPathways <- fgseaLocal[head(order(pval), n=15), pathway]
        # plotMyGseaTable("Cell Sublocation",SublocationSet2[sigPathways], ranks3, fgseaLocal, gseaParam = 0.1)### need to change the plot title "Pathway"
        plotEnrichment(SublocationSet2[["Exosomoe"]],ranks3) + labs(title="Exosome")
        fwrite(fgseaLocal, file=paste0(outF1,"Dis",diseaseGroup,filterType,colnames(sigPrTableUni)[pColum],"detailSublocationEnrich.txt"), sep="\t", sep2=c("", " ", ""))
        
        ### enrichment test for broad level of cell locations
        fgseaBroadLoc <- fgsea(pathways = SubLocationBroadSet2, stats=ranks3, scoreType = "std",eps = 0.0,minSize=10, maxSize=1000)
        # plotMyGseaTable("Cell Broad Sublocation",SubLocationBroadSet2, ranks3, fgseaBroadLoc, gseaParam = 0.5)### need to change the plot title "Pathway"
        fwrite(fgseaBroadLoc, file=paste0(outF1,"Dis",diseaseGroup,filterType,colnames(sigPrTableUni)[pColum],"BroadSubLocationEnrich.txt"), sep="\t", sep2=c("", " ", ""))
        
        ### cell/tissue type enrichment #########################################
        ### Tissue Type enrichment test
        fgseaTissue <- fgsea(pathways = TissueSet2, stats=ranks3, scoreType = "std",eps = 0.0,minSize=1, maxSize=1000)
        # plotMyGseaTable("Tissue type",TissueSet2, rank3, fgseaTissue, gseaParam = 0.2)### need to change the plot title "Pathway"
        fwrite(fgseaTissue, file=paste0(outF1,"Dis",diseaseGroup,filterType,colnames(sigPrTableUni)[pColum],"TissueEnrich.txt"), sep="\t", sep2=c("", " ", ""))
        
        ### detailed Cell Type enrichment test
        fgseaCell <- fgsea(pathways = CellSet2, stats=ranks3, scoreType = "std",eps = 0.0,minSize=10, maxSize=1000)
        # plotMyGseaTable("Cell Type",CellSet2, ranks3, fgseaCell, gseaParam=1)
        fwrite(fgseaCell, file=paste0(outF1,"Dis",diseaseGroup,filterType,colnames(sigPrTableUni)[pColum],"CellEnrich.txt"), sep="\t", sep2=c("", " ", ""))
        
        ### Enrichment test using scRNA seq data.##########
        ### perform cell type enrichment using specific cell type marker gene sets derived from OA paper
        fgseaSCRNAcellType <- fgsea(pathways = scRNAgeneCellSets2, stats=ranks3, scoreType = "std", eps = 0.0,minSize=5, maxSize=1000)
        # plotMyGseaTable("Cell Type (scRNA data)",scRNAgeneCellSets2, ranks3, fgseaSCRNAcellType, gseaParam = 0.1)### need to change the plot title "Pathway"
        fwrite(fgseaSCRNAcellType, file=paste0(outF1,"Dis",diseaseGroup,filterType,colnames(sigPrTableUni)[pColum],"scRNACellEnrich.txt"), sep="\t", sep2=c("", " ", ""))
        
        ### Enrichment of Markers of necrosis and apoptosis##########
        NecApopAlarEnrich <- fgsea(pathways = NecApopAlarGeneSet, stats=ranks1, scoreType = "std",eps = 0.0,minSize=5, maxSize=1000)
        # plotMyGseaTable("Necrosis|Apoptosis|Alarming",NecApopAlarGeneSet, ranks1, NecApopAlarEnrich, gseaParam = 0.1)
        fwrite(NecApopAlarEnrich, file=paste0(outF1,"Dis",diseaseGroup,filterType,colnames(sigPrTableUni)[pColum],"NecApopAlarEnrich.txt"), sep="\t", sep2=c("", " ", ""))
        
        ### sublocation bar plot
        ### define differentially expressed proteins for endotype
        DEtestPvalue <- p.adjust(rankPro, method = "bonferroni", n=nrow(sigPrTableUni))
        DEentrezName <- sigPrTableUni$EnrichEntrezSymbolLabel[which(DEtestPvalue<0.05&fc>1)]
        if(length(DEentrezName)==0){
          writeLines("No differentially expressed genes are found \n",log1.3)
        }else{
          
          ### bar plot for sublocation enrichment
          datFplot31 <- barComposition(fgseaLocal,SublocationSet2,1,DEentrezName,"SubLocationEnrichType","Protein.Membership.Counts","Composition of subcelluar locations of Protein signatures for endotypes")
          ### save the data frame of the enrichment plot
          write.table(datFplot31,paste0(outF1,"Dis",diseaseGroup,filterType,colnames(sigPrTableUni)[pColum],".detailSublocationComposition.csv"),row.names = FALSE,col.names=TRUE,quote=FALSE,sep="\t")
          
          ### bar plot for broad cell sublocation enrichment
          datFplot32 <- barComposition(fgseaBroadLoc,SubLocationBroadSet2,1,DEentrezName,"Broad SubLocation Enrich Type","Protein.Membership.Counts","Composition of borad subcelluar locations of Protein signatures for endotypes")
          write.table(datFplot32,paste0(outF1,"Dis",diseaseGroup,filterType,colnames(sigPrTableUni)[pColum],".SubLocationBroadComposition.csv"),row.names = FALSE,col.names=TRUE,quote=FALSE,sep="\t")
          
          ### bar plot for tissue type enrichment
          datFplot4 <- barComposition(fgseaTissue,TissueSet2,1,DEentrezName,"Tissue Type Enrich Type","Tissue.Membership.Counts","Composition of tissue types of Protein signatures for endotypes")
          write.table(datFplot4,paste0(outF1,"Dis",diseaseGroup,filterType,colnames(sigPrTableUni)[pColum],".TissueComposition.csv"),row.names = FALSE,col.names=TRUE,quote=FALSE,sep="\t")
          
          ### bar plot for cell type enrichment
          datFplot51 <- barComposition(fgseaCell,CellSet2,1,DEentrezName,"Cell Type Enrich Type","Cell.Membership.Counts","Composition of cell types of Protein signatures for endotypes")
          write.table(datFplot51,paste0(outF1,"Dis",diseaseGroup,filterType,colnames(sigPrTableUni)[pColum],".allCellComposition.csv"),row.names = FALSE,col.names=TRUE,quote=FALSE,sep="\t")
          
          ### bar plot for scRNA cell type enrichment
          datFplot52 <- barComposition(fgseaSCRNAcellType,scRNAgeneCellSets2,1,DEentrezName,"Cell Type Enrich Type","Cell.Membership.Counts","Composition of cell types of Protein signatures for endotypes (scRNA)")
          write.table(datFplot52,paste0(outF1,"Dis",diseaseGroup,filterType,colnames(sigPrTableUni)[pColum],".scRNAgeneCellComposition.csv"),row.names = FALSE,col.names=TRUE,quote=FALSE,sep="\t")
          
          
          ###GWAS endo heritatiblity enrichment test for OA riks variants###################
          ### get ensemble genes for LDSC analysis in R.
          ### inside pCol loop, generate endo gene set.Cntrol gene set has been generated before pCol loop begins
          DEentrezName2 <- sigPrTableUni$EntrezGeneSymbol[which(DEtestPvalue<0.05)] ### the following need original entrez infomation before protein complex replacement
          HUGOs <- unique(getBM(attributes=c("hgnc_symbol"),filters = "entrezgene_accession", values=DEentrezName2, mart= humGenes)[,]) ###hgnc used to retrevie chromo position
          
          DEensemble <- getBM(attributes=c("ensembl_gene_id","chromosome_name","start_position", "end_position"),filters = "hgnc_symbol", values=HUGOs, mart= humGenes)
          DEensemble2 <- DEensemble[which(!grepl("_|X",DEensemble$chromosome_name)),]
          colnames(DEensemble2) = c("GENE","CHR","START","END") ### consistent with ldsc python "gene_set_to_bed" function
          DEensemble2$CHR <- paste0("ch",DEensemble2$CHR)
          DEensembleName <- as.matrix(DEensemble2[,1],ncol=1)
          
          write.table(DEensembleName,paste0(outF2,"/DEPcol.",pColum,".GeneSet"),row.names = FALSE,col.names=FALSE,quote=FALSE)
          write.table(DEensemble2,paste0(outF2,"/DEPcol.",pColum,".Location.txt"),row.names = FALSE,col.names=TRUE,quote=FALSE)
          writeLines(paste0("Onput file:",paste0(outF2,"/DEPcol.",pColum,".GeneSet"),"\n"),log1.3)
          writeLines(paste0("Onput file:",paste0(outF2,"/DEPcol.",pColum,".Location.txt"),"\n"),log1.3)
        }
        # ### to run LDSC, not nessesary to install pacakges individually. Install Anaconda. Need "conda env create --file environment.yml" once, but need "source activate ldsc" each time run the package.
        # ### panda download via pip. After Anoconda is installed, "conda env create --file environment.yml" report "Cython-generated file 'pandas/_libs/testing.c' not found.Cython is required to compile pandas from a development branch.Please install Cython or download a release package of pandas.", so we do "curl https://bootstrap.pypa.io/get-pip.py | python3" "pip install pandas"
        # ### run bash when get the ldsc working envinronment ready. bash l2for22.sh. l2 calculation loop through 22 chromosomes + heritability enrichment test
        # ### output endo1.results.txt
        # 
        ### upstream analysis
        ### sigPrTableUni is for unique somaID, upstream analysis, we need unique ENTREZ id.
        sigPrTableUni2 = sigPrTableUni[sapply(unique(sigPrTableUni[,"EntrezGeneID"]),function(x){which(sigPrTableUni[,"EntrezGeneID"]==x)[1]}),]
        EntrezCol = which(colnames(sigPrTableUni2)=="EntrezGeneID")
        fcCol = which(colnames(sigPrTableUni2)=="fc")
        upstreamDat = sigPrTableUni2[,c(EntrezCol,fcCol,pColum)] ### columns of sigProTable is not constant, so need to dynamically find the relevant columns
        colnames(upstreamDat) <- c("entrez", "fc", "pvalue") ### name required by the package, package require unique gene ID
        upstreamDat$fc = as.numeric(upstreamDat$fc)
        upstreamDat$pvalue = as.numeric(upstreamDat$pvalue)

        try(quaternary_results <- QuaternaryProd::RunCRE_HSAStringDB(upstreamDat[which(upstreamDat$entrez!=''),], method = "Quaternary",
                                                                     fc.thresh = log2(1.3), pval.thresh = 0.05/nrow(upstreamDat),
                                                                     only.significant.pvalues = TRUE,
                                                                     significance.level = 0.05,
                                                                     epsilon = 1e-16, progressBar = FALSE,
                                                                     relations = NULL, entities = NULL))

        if(exists("quaternary_results")){
          goodPid = which(quaternary_results[,"pvalue"]!=-1 & quaternary_results[,"pvalue"]<0.05/nrow(upstreamDat) & quaternary_results[,"symbol"]!="No-Symbol")
          print(paste(length(goodPid),"significant upstream regulators are found."))

          ###Top 10 regulators
          topRegulators = quaternary_results[goodPid, c("uid","symbol","regulation","pvalue")]
          x= list("UP Regulation" = which(quaternary_results[goodPid,"regulation"]=="up"),"Down Regulation" = which(quaternary_results[goodPid,"regulation"]=="down"))
          # ggvenn(x,fill_color = c("#0073C2FF", "#CD534CFF"),stroke_size = 0.5, set_name_size = 4)
          rm("quaternary_results")
          write.table(topRegulators,file=paste0(outF1,"Dis",diseaseGroup,filterType,".Pcol.",colnames(sigPrTableUni)[pColum],".TopRegulators.txt"))
          writeLines(paste0("Onput file:",paste0(outF1,"Dis",diseaseGroup,filterType,".Pcol.",colnames(sigPrTableUni)[pColum],"TopRegulators.txt"),"\n"),log1.3)
        }else{writeLines("no significant upstream regulators are found.\n",paste0(outF1,"Dis",diseaseGroup,filterType,"Pcol",pColum,"TopRegulators.txt"))}

        dev.off()
      }
      writeLines(paste0("p column ",pColum," fininished","\n"),log1.3)
    }
    
    ### for ldsc gwas analysis, automatically generate a "OAendo.ldcts" file for heritatilibty enrichment test
    ### when the files are not directly under ldsc folder, the whole file path should be included to initialise the parameters
    allDEsetfiles <- list.files(paste0(myPathIn5,"ldsc/DiseaseGroup",diseaseGroup,"/",filterType,"/"),"DEPcol.*.GeneSet")
    allCategories <- sub("GeneSet","",allDEsetfiles)
    allDEsetfiles2 <- paste0(myPathIn5,"ldsc/DiseaseGroup",diseaseGroup,"/",filterType,"/endoEnrich_ldscores/",allCategories,",",myPathIn5,"ldsc/DiseaseGroup",diseaseGroup,"/",filterType,"/endoEnrich_ldscores/Control.")
    
    OActs <- cbind(allCategories,allDEsetfiles2)
    write.table(OActs,paste0(myPathIn5,"ldsc/DiseaseGroup",diseaseGroup,"/",filterType,"/OAendo.ldcts"),sep="\t",row.names = FALSE, col.names=FALSE, quote=FALSE)
    
    writeLines(paste0(filterType," fininished","\n"),log1.3)
    
  }
  writeLines(paste0("Disease group ",diseaseGroup," fininished","\n"),log1.3)
}

writeLines(paste0("finish time ",Sys.Date(),".\n"),log1.3)
#if run in R studio, sink the warnings to log1.3. 
sink(paste0(myPathOut,"1.3.log"), open = "wt", append = TRUE, type = c("message"),split = FALSE)
sink()
close(log1.3)

### reset output to log file back to console window for future usage.
closeAllConnections()
#in bash: add this: Rscript /Users/ydeng/Documents//DA__202111/DAcode2021.11/PrimaryAnalysis/subanalysis1.3.R 2> RWarning.txt
