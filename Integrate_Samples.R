## to replicate conda environment used in paper, install conda, and then call the following in terminal: conda env create -f environment_droplet.yml - yml file here: https://github.com/brianherb/HumanHypothalamusDev/blob/main/environment_droplet.yml

library(Seurat)
library(limma)

devtools::source_url('https://github.com/brianherb/HumanHypothalamusDev/blob/main/HypoPub_functions_reference.R?raw=TRUE')

######################
### load datasets ####
######################

## Lein hypo samples 
load(url('https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/EdHypoSamples_PostSCT.rda')) #EdHypoNormDat

## Kriegstein data 
load(url('https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/HypoSamples_PostSCT.rda')) #HypoNormDat

## Non-hypothalamus Human Embryonic samples 

## GW18 Cortex
load(url('https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/GW18CtxSamples_PostSCT.rda')) #GW18CtxNormDat

## GW18 GE
load(url('https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/GW18GESamples_PostSCT.rda')) #GW18GENormDat

## GW19
load(url('https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/GW19_NonHySamples_PostSCT.rda')) # GW19_NormDat

## GW20
load(url('https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/GW20_34_NonHySamples_PostSCT.rda')) # GW20_34_NormDat


############################################
### Integrate CS13 - GW25, all cells    ####
############################################

combDat = HypoNormDat[c("CS13","CS14","CS15","CS22_hypo","CS22_2_hypo","GW16_hypo","GW18_hypo","GW19_hypo","GW20_34_hypo","GW22T_hypo1","GW25_3V_hypo")]  
IntName = 'HypoCS13_GW25_mt10'  
PrintColCat = c('sample','seurat_clusters') 
PrintColCont = c('nCount_RNA','nFeature_RNA','percent.mt','GAD1','GAD2','NEUROD6','SLC17A6')
#40927 cells in integrated object

## Integrated object saved here: 'https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/HypoCS13_GW25_mt10_integrated.rds'

######################################################
### Integrate All Lein adult samples, all cells   ####
######################################################

combDat = EdHypoNormDat  
IntName = 'EdHypoAll_mt10'  
PrintColCat = c('sample','seurat_clusters') 
PrintColCont = c('nCount_RNA','nFeature_RNA','percent.mt','GAD1','GAD2','SLC17A6')
#35338 cells in integrated object

## Integrated object saved here: 'https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/EdHypoAll_mt10_integrated.rds'

######################################################
### Integrate all Lein adult samples, neurons     ####
######################################################

## this is done after creating 'EdHypoAll_mt10' - using surat groups that express GAD1/2 and SLC17A6

combDat = EdHypoNormDat  

exInd = unique(which(is.na(match(EdHypoAll_mt10@meta.data$seurat_clusters,c(2:4,6,9,28,30:33))))) 

EdNeuronBarcodes_mt10 = colnames(EdHypoAll_mt10)[exInd]

save(EdNeuronBarcodes_mt10,file='./SeuratObj/EdNeuronBarcodes_mt10.rda',compress=TRUE) #25424

## EdNeuronBarcodes_mt10 saved here: 'https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/EdNeuronBarcodes_mt10.rda'

for(i in names(combDat)){
combDat[[i]] = subset(combDat[[i]],cells=intersect(colnames(combDat[[i]]),EdNeuronBarcodes_mt10))
cat(paste(i,', ',sep=''))
}

IntName = 'EdHypoNeuro_mt10'  
PrintColCat = c('sample') 
PrintColCont = c('nCount_RNA','nFeature_RNA','percent.mt','GAD1','GAD2','SLC17A6')

## Integrated object saved here: 'https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/EdHypoNeuro_mt10_integrated.rds'


###########################################################
### Integrate all Kriegstein neurons, starting at CS22 ####
###########################################################

## this is done after creating 'HypoCS13_GW25_mt10' - using seurat groups that express GAD1/2 and SLC17A6 (and are more mature - ignoring progenitors )

combDat = HypoNormDat[c("CS22_hypo","CS22_2_hypo","GW16_hypo","GW18_hypo","GW19_hypo","GW20_34_hypo","GW22T_hypo1","GW25_3V_hypo")] 

KaNeuronBarcodes_mt10 = colnames(HypoCS13_GW25_mt10)[which(!is.na(match(HypoCS13_GW25_mt10@meta.data$seurat_clusters,c(0,4,15,19,20,25))))]

save(KaNeuronBarcodes_mt10,file='./SeuratObj/KaNeuronBarcodes_mt10.rda',compress=TRUE)

## KaNeuronBarcodes_mt10 saved here: 'https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/KaNeuronBarcodes_mt10.rda'

for(i in names(combDat)){
combDat[[i]] = subset(combDat[[i]],cells=intersect(colnames(combDat[[i]]),KaNeuronBarcodes_mt10))
cat(paste(i,', ',sep=''))
}

IntName = 'KaHypoNeuro_mt10'  
PrintColCat = c('sample','seurat_clusters') 
PrintColCont = c('nCount_RNA','nFeature_RNA','percent.mt','GAD1','GAD2','SLC17A6')

## Integrated object saved here: 'https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/KaHypoNeuro_mt10_integrated.rds'

##############################################
### Integrate Kriegstein and Lein neurons ####
##############################################

combDatKA = HypoNormDat[c("CS22_hypo","CS22_2_hypo","GW16_hypo","GW18_hypo","GW19_hypo","GW20_34_hypo","GW22T_hypo1","GW25_3V_hypo")] 

for(i in names(combDatKA)){
combDatKA[[i]] = subset(combDatKA[[i]],cells=intersect(colnames(combDatKA[[i]]),KaNeuronBarcodes_mt10))
cat(paste(i,', ',sep=''))
}

combDatED = EdHypoNormDat  

for(i in names(combDatED)){
combDatED[[i]] = subset(combDatED[[i]],cells=intersect(colnames(combDatED[[i]]),EdNeuronBarcodes_mt10))
cat(paste(i,', ',sep=''))
}

combDat = c(combDatKA,combDatED)

IntName = 'EdKaHypoNeurons_mt10'  
PrintColCat = c('sample','seurat_clusters') 
PrintColCont = c('nCount_RNA','nFeature_RNA','percent.mt','GAD1','GAD2','SLC17A6')

## Integrated object saved here: 'https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/EdKaHypoNeurons_mt10_3D_integrated.rds'


##############################################
### Integrate Ctx GE Hypo 						####
##############################################

## Done in Fig4_HumanHypoCtxGE.R 

############################################
### Common code for integration 	    ####
############################################

for (i in 1:length(combDat)) {
	DefaultAssay(combDat[[i]]) = 'RNA'
	combDat[[i]] <- subset(combDat[[i]], percent.mt<=10)
	if('FOXG1' %in% rownames(combDat[[i]]@assays$RNA) ){
		combDat[[i]] <- subset(combDat[[i]], cells=colnames(combDat[[i]])[which(combDat[[i]]@assays$RNA@counts['FOXG1',]==0)])
	}
	if('NEUROD6' %in% rownames(combDat[[i]]@assays$RNA) ){
		combDat[[i]] <- subset(combDat[[i]], cells=colnames(combDat[[i]])[which(combDat[[i]]@assays$RNA@counts['NEUROD6',]==0)])
	}
    combDat[[i]] <- NormalizeData(combDat[[i]], verbose = FALSE)
    combDat[[i]] <- FindVariableFeatures(combDat[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    cat(paste(i,', ',sep=''))
}

features <- SelectIntegrationFeatures(object.list = combDat)

for (i in 1:length(combDat)) {
 combDat[[i]]<- ScaleData(combDat[[i]], features = features, verbose = FALSE)
    combDat[[i]] <- RunPCA(combDat[[i]], features = features, verbose = FALSE)
}

All.anchors <- FindIntegrationAnchors(object.list = combDat, anchor.features = features, reduction = "rpca") 
All.integrated <- IntegrateData(anchorset = All.anchors, dims = 1:30) # was k.weight =30 dims =1:15 
DefaultAssay(All.integrated) <- "integrated"
All.integrated <- ScaleData(All.integrated, verbose = FALSE)
All.integrated <- RunPCA(All.integrated, npcs = 30, verbose = FALSE)
All.integrated <- RunUMAP(All.integrated, reduction = "pca", dims = 1:30)

All.integrated <- FindNeighbors(All.integrated, dims = 1:30, verbose = FALSE)
All.integrated <- FindClusters(All.integrated, verbose = FALSE, resolution=0.5)
saveRDS(All.integrated,file=paste('./SeuratObj/',IntName,'_integrated.rds',sep=''),compress=TRUE) 

DefaultAssay(All.integrated) = 'RNA'

pdf(file=paste("./TestPlots/",IntName,"_Check_Clustering.pdf",sep=''),width=12,height=8)
for(k in PrintColCat){
print(DimPlot(All.integrated, reduction = "umap", group.by = k, label = TRUE, repel = TRUE))
}
for(m in PrintColCont){
print(FeaturePlot(All.integrated, features = m, pt.size = 0.2))
}

dev.off()