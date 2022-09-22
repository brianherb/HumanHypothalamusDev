## to replicate conda environment used in paper, install conda, and then call the following in terminal: conda env create -f environment_droplet.yml - yml file here: https://github.com/brianherb/HumanHypothalamusDev/blob/main/environment_droplet.yml

library(SingleCellExperiment)
library(Seurat)
library(scDblFinder)
library(devtools)

devtools::source_url('https://github.com/brianherb/HumanHypothalamusDev/blob/main/HypoPub_functions_reference.R?raw=TRUE')

## Paths in these scripts mirror directories used in https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/ We suggest that you establish a single working directory and add the following folders: Analysis, Counts, Networks, Reference_Files, SeuratObj, and TestPlots - a directory can be created in R like so: mkdir('./Counts')

setwd('./Counts')

## download counts from https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/Counts/

system('wget -r -np -nH --cut-dirs=12 -e robots=off --reject="index.html*" https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/Counts/')
## alternatively, you can download the data from NeMO Archive - https://nemoarchive.org

setwd('./')


HypoSamples = c("CS13","CS14","CS15","CS22_hypo","CS22_2_hypo","GW16_hypo","GW18_hypo","GW19_hypo","GW20_34_hypo","GW22T_hypo1","GW25_3V_hypo")

CS13_data <- Read10X(data.dir = "./Counts/CS13prosencephalon_GRCh38")
CS14_data <- Read10X(data.dir = "./Counts/CS14cortex_GRCh38")
CS15_data <- Read10X(data.dir = "./Counts/CS15_forebrain") 
CS22_hypo_data <- Read10X(data.dir = "./Counts/CS22_Hypothalamus" )
CS22_2_hypo_data <- Read10X(data.dir = "./Counts/CS22_2_hypothalamus")
GW16_hypo_data <- Read10X(data.dir =  "./Counts/GW16_hypo")      
GW18_hypo_data <- Read10X(data.dir = "./Counts/GW18_hypothalamus") 
GW19_hypo_data <- Read10X(data.dir =  "./Counts/GW19_hypothalamus") 
GW20_34_hypo_data <- Read10X(data.dir = "./Counts/GW20_34_hypo")     
GW22T_hypo1_data <- Read10X(data.dir =  "./Counts/GW22T_hypo1")  
GW25_3V_hypo_data <- Read10X(data.dir =  "./Counts/GW25_3V_hypo")


HypoRawDat = list(CS13_data,CS14_data,CS15_data,CS22_hypo_data,CS22_2_hypo_data,GW16_hypo_data,GW18_hypo_data,GW19_hypo_data,GW20_34_hypo_data,GW22T_hypo1_data,GW25_3V_hypo_data)

names(HypoRawDat) = HypoSamples

HypoCellTracking = data.frame(sample=HypoSamples,Raw=0,AfterCutoffs=0,AfterDoublets=0)
rownames(HypoCellTracking) = HypoSamples

HypoDat = vector(mode='list',length=length(HypoSamples))
names(HypoDat) = HypoSamples
HypoNormDat = vector(mode='list',length=length(HypoSamples))
names(HypoNormDat) = HypoSamples


for(i in HypoSamples){

HypoDat[[i]] <- CreateSeuratObject(counts = HypoRawDat[[i]])
HypoCellTracking[i,'Raw'] = ncol(HypoDat[[i]])
HypoDat[[i]]@meta.data$sample=i
HypoDat[[i]] = RenameCells(HypoDat[[i]],add.cell.id=i)
HypoDat[[i]] <- PercentageFeatureSet(HypoDat[[i]], pattern = "^MT-", col.name = "percent.mt")
## QC plots
pdf(file=paste("./TestPlots/",i,"_all_QC.pdf",sep=''),width=12,height=8)
print(VlnPlot(HypoDat[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
)
print(FeatureScatter(HypoDat[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt"))
print(FeatureScatter(HypoDat[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))
dev.off()
HypoNormDat[[i]] = subset(HypoDat[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA > 1000& nCount_RNA < 15000 & percent.mt < 20)
HypoCellTracking[i,'AfterCutoffs'] = ncol(HypoNormDat[[i]])
if(ncol(HypoNormDat[[i]])<100) next
HypoNormDat[[i]] <- SCTransform(HypoNormDat[[i]], vars.to.regress = "percent.mt", verbose = FALSE)
HypoNormDat[[i]] <- RunPCA(HypoNormDat[[i]], verbose = FALSE)
HypoNormDat[[i]] <- RunUMAP(HypoNormDat[[i]], dims = 1:50, verbose = FALSE)
HypoNormDat[[i]] <- FindNeighbors(HypoNormDat[[i]], dims = 1:50, verbose = FALSE)
HypoNormDat[[i]] <- FindClusters(HypoNormDat[[i]], verbose = FALSE, resolution=0.5)

DefaultAssay(HypoNormDat[[i]]) = 'RNA'

pdf(file=paste("./TestPlots/",i,"_Check_Cluster_PreDupRm.pdf",sep=''),width=12,height=8)
print(DimPlot(HypoNormDat[[i]], reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE))
print(FeaturePlot(HypoNormDat[[i]], features = c("AQP4","AGT"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(HypoNormDat[[i]], features = c("FREM2","FZD2"), pt.size = 0.2, ncol = 2) )
print(FeaturePlot(HypoNormDat[[i]], features = c("ITM2A","IGFBP7"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(HypoNormDat[[i]], features = c("MATN4","SCRG1"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(HypoNormDat[[i]], features = c("SMAGP","HPGD"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(HypoNormDat[[i]], features = c("ERMN","MOBP"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(HypoNormDat[[i]], features = c("SLC44A1","GPR17"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(HypoNormDat[[i]], features = c("C1QB","HEXB"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(HypoNormDat[[i]], features = c("SLC17A6","GAD1"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(HypoNormDat[[i]], features = c("AQP4","RAX"), pt.size = 0.2, ncol = 2))
dev.off()

## doublet id and removal -  
DefaultAssay(HypoNormDat[[i]]) = 'SCT'
Seurat_Object_Diet <- DietSeurat(HypoNormDat[[i]], graphs = "pca")
TMPsce <- as.SingleCellExperiment(Seurat_Object_Diet)
#TMPsce = as.SingleCellExperiment(HypoNormDat[[i]])
TMPsce <- scDblFinder(TMPsce)
HypoNormDat[[i]]@meta.data$scDblFinder = colData(TMPsce)$scDblFinder.class
HypoNormDat[[i]] = subset(HypoDat[[i]],cells = colnames(HypoNormDat[[i]])[which(HypoNormDat[[i]]@meta.data$scDblFinder=='singlet')]) #7894
HypoCellTracking[i,'AfterDoublets'] = ncol(HypoNormDat[[i]])
if(ncol(HypoNormDat[[i]])<100) next
HypoNormDat[[i]] <- SCTransform(HypoNormDat[[i]], vars.to.regress = "percent.mt", verbose = FALSE)
HypoNormDat[[i]] <- RunPCA(HypoNormDat[[i]], verbose = FALSE)
HypoNormDat[[i]] <- RunUMAP(HypoNormDat[[i]], dims = 1:50, verbose = FALSE)
HypoNormDat[[i]] <- FindNeighbors(HypoNormDat[[i]], dims = 1:50, verbose = FALSE)
HypoNormDat[[i]] <- FindClusters(HypoNormDat[[i]], verbose = FALSE, resolution=0.5)

pdf(file=paste("./TestPlots/",i,"_Check_Cluster_PostDupRm.pdf",sep=''),width=12,height=8)
print(DimPlot(HypoNormDat[[i]], reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE))
print(FeaturePlot(HypoNormDat[[i]], features = c("AQP4","AGT"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(HypoNormDat[[i]], features = c("FREM2","FZD2"), pt.size = 0.2, ncol = 2) )
print(FeaturePlot(HypoNormDat[[i]], features = c("ITM2A","IGFBP7"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(HypoNormDat[[i]], features = c("MATN4","SCRG1"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(HypoNormDat[[i]], features = c("SMAGP","HPGD"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(HypoNormDat[[i]], features = c("ERMN","MOBP"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(HypoNormDat[[i]], features = c("SLC44A1","GPR17"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(HypoNormDat[[i]], features = c("C1QB","HEXB"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(HypoNormDat[[i]], features = c("SLC17A6","GAD1"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(HypoNormDat[[i]], features = c("AQP4","RAX"), pt.size = 0.2, ncol = 2))
dev.off()

cat(paste('\n\n',i,' done normalization','\n\n',sep=''))

}

save(HypoDat,file='./SeuratObj/HypoSamples_AllCells.rda',compress=TRUE)

## used in paper:
## load(url('https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/HypoSamples_AllCells.rda'))

save(HypoNormDat,file='./SeuratObj/HypoSamples_PostSCT.rda',compress=TRUE)

## used in paper:
## load(url('https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/HypoSamples_PostSCT.rda'))

write.csv(HypoCellTracking,file='./SeuratObj/HypoSamples_CellCounts.csv')


## Cortex and GE samples from GW18, GW19 and GW20

## GW18 

cortexMeta = read.csv('https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/Reference_Files/Primary_samples_cortex.csv') ## Cell type assignments from the work leading towards the publication "An atlas of cortical arealization identifies dynamic molecular signatures" Bhaduri 2021

cortexCells = flexsplit(cortexMeta$Cell,'_')[,1]

CortexSamples = c("GW18_V1","GW18_parietal","GW18_motor","GW18_PFC","GW18_somato","GW18_somatosensory","GW18_temporal")

GW18_V1_data <- Read10X(data.dir = "./Counts/GW18_V1") 
GW18_parietal_data <- Read10X(data.dir = "./Counts/GW18_parietal") 
GW18_motor_data <- Read10X(data.dir = "./Counts/GW18_motor") 
GW18_PFC_data <- Read10X(data.dir = "./Counts/GW18_PFC") 
GW18_somato_data <- Read10X(data.dir = "./Counts/GW18_somato")

length(intersect(cortexCells,colnames(GW18_V1_data))) #10758
length(intersect(cortexCells,colnames(GW18_parietal_data))) #8786
length(intersect(cortexCells,colnames(GW18_motor_data))) #16954
length(intersect(cortexCells,colnames(GW18_PFC_data))) #15108
length(intersect(cortexCells,colnames(GW18_somato_data))) # 11080

GW18CtxMeta=cortexMeta[which(cortexMeta$Individual=='GW18'),]

GW18CtxMeta$sample = paste(GW18CtxMeta$Individual,GW18CtxMeta$Area,sep='_')

GW18CtxSamples = c("GW18_V1","GW18_parietal","GW18_motor","GW18_PFC","GW18_somatosensory")

GW18CtxRawDat = list(GW18_V1_data,GW18_parietal_data,GW18_motor_data,GW18_PFC_data,GW18_somato_data)
names(GW18CtxRawDat) = GW18CtxSamples

## cell number tracking 

GW18CtxCellTracking = data.frame(sample=GW18CtxSamples,Raw=0,AfterCutoffs=0,AfterDoublets=0)
rownames(GW18CtxCellTracking) = GW18CtxSamples

GW18CtxAllDat = vector(mode='list',length=length(GW18CtxSamples))
names(GW18CtxAllDat) = GW18CtxSamples
GW18CtxNormDat = vector(mode='list',length=length(GW18CtxSamples))
names(GW18CtxNormDat) = GW18CtxSamples


for(i in GW18CtxSamples){

GW18CtxAllDat[[i]] <- CreateSeuratObject(counts = GW18CtxRawDat[[i]])
GW18CtxCellTracking[i,'Raw'] = ncol(GW18CtxAllDat[[i]])
GW18CtxAllDat[[i]]@meta.data$sample=i
GW18CtxAllDat[[i]] = RenameCells(GW18CtxAllDat[[i]],add.cell.id=i)
GW18CtxAllDat[[i]] <- PercentageFeatureSet(GW18CtxAllDat[[i]], pattern = "^MT-", col.name = "percent.mt")
## QC plots
pdf(file=paste("./TestPlots/",i,"_all_QC.pdf",sep=''),width=12,height=8)
print(VlnPlot(GW18CtxAllDat[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
)
print(FeatureScatter(GW18CtxAllDat[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt"))
print(FeatureScatter(GW18CtxAllDat[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))
dev.off()

tmpMeta = GW18CtxMeta[which(GW18CtxMeta$sample==i),]
tmpMetaCell=flexsplit(tmpMeta$Cell,'_')[,1]
tmpMetaCell = paste(i,tmpMetaCell,sep='_')
rownames(tmpMeta) = tmpMetaCell
GW18CtxAllDat[[i]]@meta.data$Cell_Class = tmpMeta[colnames(GW18CtxAllDat[[i]]),'Class']
GW18CtxAllDat[[i]]@meta.data$Cell_State = tmpMeta[colnames(GW18CtxAllDat[[i]]),'State']
GW18CtxAllDat[[i]]@meta.data$Cell_Type = tmpMeta[colnames(GW18CtxAllDat[[i]]),'Type']

GW18CtxNormDat[[i]] = subset(GW18CtxAllDat[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA > 1000& nCount_RNA < 15000 & percent.mt < 20)
GW18CtxCellTracking[i,'AfterCutoffs'] = ncol(GW18CtxNormDat[[i]])
GW18CtxNormDat[[i]] <- SCTransform(GW18CtxNormDat[[i]], vars.to.regress = "percent.mt", verbose = FALSE)
GW18CtxNormDat[[i]] <- RunPCA(GW18CtxNormDat[[i]], verbose = FALSE)
GW18CtxNormDat[[i]] <- RunUMAP(GW18CtxNormDat[[i]], dims = 1:50, verbose = FALSE)
GW18CtxNormDat[[i]] <- FindNeighbors(GW18CtxNormDat[[i]], dims = 1:50, verbose = FALSE)
GW18CtxNormDat[[i]] <- FindClusters(GW18CtxNormDat[[i]], verbose = FALSE, resolution=0.5)

DefaultAssay(GW18CtxNormDat[[i]]) = 'RNA'

pdf(file=paste("./TestPlots/",i,"_Check_Cluster_PreDupRm.pdf",sep=''),width=12,height=8)
DimPlot(GW18CtxNormDat[[i]], label = TRUE) + NoLegend()
DimPlot(GW18CtxNormDat[[i]], reduction = "umap", group.by = "Cell_Class", label = TRUE, repel = TRUE)
DimPlot(GW18CtxNormDat[[i]], reduction = "umap", group.by = "Cell_State", label = TRUE, repel = TRUE)
DimPlot(GW18CtxNormDat[[i]], reduction = "umap", group.by = "Cell_Type", label = TRUE, repel = TRUE)
print(FeaturePlot(GW18CtxNormDat[[i]], features = c("AQP4","AGT"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW18CtxNormDat[[i]], features = c("FREM2","FZD2"), pt.size = 0.2, ncol = 2) )
print(FeaturePlot(GW18CtxNormDat[[i]], features = c("ITM2A","IGFBP7"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW18CtxNormDat[[i]], features = c("MATN4","SCRG1"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW18CtxNormDat[[i]], features = c("SMAGP","HPGD"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW18CtxNormDat[[i]], features = c("ERMN","MOBP"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW18CtxNormDat[[i]], features = c("SLC44A1","GPR17"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW18CtxNormDat[[i]], features = c("C1QB","HEXB"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW18CtxNormDat[[i]], features = c("SLC17A6","GAD1"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW18CtxNormDat[[i]], features = c("AQP4","RAX"), pt.size = 0.2, ncol = 2))
dev.off()

## doublet id and removal -  
DefaultAssay(GW18CtxNormDat[[i]]) = 'SCT'
TMPsce = Seurat::as.SingleCellExperiment(GW18CtxNormDat[[i]])
TMPsce <- scDblFinder(TMPsce,clust.method='overcluster')
GW18CtxNormDat[[i]]@meta.data$scDblFinder = colData(TMPsce)$scDblFinder.class
GW18CtxNormDat[[i]] = subset(GW18CtxAllDat[[i]],cells = colnames(GW18CtxNormDat[[i]])[which(GW18CtxNormDat[[i]]@meta.data$scDblFinder=='singlet')]) #7894
GW18CtxCellTracking[i,'AfterDoublets'] = ncol(GW18CtxNormDat[[i]])
GW18CtxNormDat[[i]] <- SCTransform(GW18CtxNormDat[[i]], vars.to.regress = "percent.mt", verbose = FALSE)
GW18CtxNormDat[[i]] <- RunPCA(GW18CtxNormDat[[i]], verbose = FALSE)
GW18CtxNormDat[[i]] <- RunUMAP(GW18CtxNormDat[[i]], dims = 1:50, verbose = FALSE)
GW18CtxNormDat[[i]] <- FindNeighbors(GW18CtxNormDat[[i]], dims = 1:50, verbose = FALSE)
GW18CtxNormDat[[i]] <- FindClusters(GW18CtxNormDat[[i]], verbose = FALSE, resolution=0.5)

pdf(file=paste("./TestPlots/",i,"_Check_Cluster_PostDupRm.pdf",sep=''),width=12,height=8)
print(DimPlot(GW18CtxNormDat[[i]], label = TRUE) + NoLegend())
print(DimPlot(GW18CtxNormDat[[i]], reduction = "umap", group.by = "Cell_Class", label = TRUE, repel = TRUE))
print(DimPlot(GW18CtxNormDat[[i]], reduction = "umap", group.by = "Cell_State", label = TRUE, repel = TRUE))
print(DimPlot(GW18CtxNormDat[[i]], reduction = "umap", group.by = "Cell_Type", label = TRUE, repel = TRUE))
print(FeaturePlot(GW18CtxNormDat[[i]], features = c("AQP4","AGT"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW18CtxNormDat[[i]], features = c("FREM2","FZD2"), pt.size = 0.2, ncol = 2) )
print(FeaturePlot(GW18CtxNormDat[[i]], features = c("ITM2A","IGFBP7"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW18CtxNormDat[[i]], features = c("MATN4","SCRG1"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW18CtxNormDat[[i]], features = c("SMAGP","HPGD"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW18CtxNormDat[[i]], features = c("ERMN","MOBP"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW18CtxNormDat[[i]], features = c("SLC44A1","GPR17"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW18CtxNormDat[[i]], features = c("C1QB","HEXB"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW18CtxNormDat[[i]], features = c("SLC17A6","GAD1"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW18CtxNormDat[[i]], features = c("AQP4","RAX"), pt.size = 0.2, ncol = 2))
dev.off()

cat(paste('\n\n',i,' done normalization','\n\n',sep=''))

}

save(GW18CtxAllDat,file='./SeuratObj/GW18CtxSamples_AllCells.rda',compress=TRUE)

## used in paper:
## load(url('https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/GW18CtxSamples_AllCells.rda'))

save(GW18CtxNormDat,file='./SeuratObj/GW18CtxSamples_PostSCT.rda',compress=TRUE)

## used in paper:
## load(url('https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/GW18CtxSamples_PostSCT.rda'))

write.csv(GW18CtxCellTracking,file='./SeuratObj/GW18CtxSamples_CellCounts.csv')



pdf(file=paste("./TestPlots/all_Cortex_regions_Check_Cluster_PostDupRm.pdf",sep=''),width=12,height=8)
for(i in GW18CtxSamples){
	print(i)
print(DimPlot(GW18CtxNormDat[[i]], label = TRUE) + NoLegend())
print(DimPlot(GW18CtxNormDat[[i]], reduction = "umap", group.by = "Cell_Class", label = TRUE, repel = TRUE))
print(DimPlot(GW18CtxNormDat[[i]], reduction = "umap", group.by = "Cell_State", label = TRUE, repel = TRUE))
print(DimPlot(GW18CtxNormDat[[i]], reduction = "umap", group.by = "Cell_Type", label = TRUE, repel = TRUE))
}
dev.off()


### All GW18 GE samples 
GW18GESamples = c("GW18_CGE","GW18_LGE","GW18_MGE")

GW18_mge_data <- Read10X(data.dir = "./Counts/GW18_MGE") 
GW18_cge_data <- Read10X(data.dir = "./Counts/GW18_CGE")
GW18_lge_data <- Read10X(data.dir = "./Counts/GW18_LGE")

GW18GERawDat = list(GW18_cge_data,GW18_lge_data,GW18_mge_data)
names(GW18GERawDat) = GW18GESamples

## cell number tracking 

GW18GECellTracking = data.frame(sample=GW18GESamples,Raw=0,AfterCutoffs=0,AfterDoublets=0)
rownames(GW18GECellTracking) = GW18GESamples

GW18GEAllDat = vector(mode='list',length=length(GW18GESamples))
names(GW18GEAllDat) = GW18GESamples
GW18GENormDat = vector(mode='list',length=length(GW18GESamples))
names(GW18GENormDat) = GW18GESamples


for(i in GW18GESamples){

GW18GEAllDat[[i]] <- CreateSeuratObject(counts = GW18GERawDat[[i]])
GW18GECellTracking[i,'Raw'] = ncol(GW18GEAllDat[[i]])
GW18GEAllDat[[i]]@meta.data$sample=i
GW18GEAllDat[[i]] = RenameCells(GW18GEAllDat[[i]],add.cell.id=i)
GW18GEAllDat[[i]] <- PercentageFeatureSet(GW18GEAllDat[[i]], pattern = "^MT-", col.name = "percent.mt")
## QC plots
pdf(file=paste("./TestPlots/",i,"_all_QC.pdf",sep=''),width=12,height=8)
print(VlnPlot(GW18GEAllDat[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
)
print(FeatureScatter(GW18GEAllDat[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt"))
print(FeatureScatter(GW18GEAllDat[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))
dev.off()

GW18GENormDat[[i]] = subset(GW18GEAllDat[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA > 1000& nCount_RNA < 15000 & percent.mt < 20)
GW18GECellTracking[i,'AfterCutoffs'] = ncol(GW18GENormDat[[i]])
GW18GENormDat[[i]] <- SCTransform(GW18GENormDat[[i]], vars.to.regress = "percent.mt", verbose = FALSE)
GW18GENormDat[[i]] <- RunPCA(GW18GENormDat[[i]], verbose = FALSE)
GW18GENormDat[[i]] <- RunUMAP(GW18GENormDat[[i]], dims = 1:50, verbose = FALSE)
GW18GENormDat[[i]] <- FindNeighbors(GW18GENormDat[[i]], dims = 1:50, verbose = FALSE)
GW18GENormDat[[i]] <- FindClusters(GW18GENormDat[[i]], verbose = FALSE, resolution=0.5)

DefaultAssay(GW18GENormDat[[i]]) = 'RNA'

pdf(file=paste("./TestPlots/",i,"_Check_Cluster_PreDupRm.pdf",sep=''),width=12,height=8)
DimPlot(GW18GENormDat[[i]], label = TRUE) + NoLegend()
print(FeaturePlot(GW18GENormDat[[i]], features = c("AQP4","AGT"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW18GENormDat[[i]], features = c("FREM2","FZD2"), pt.size = 0.2, ncol = 2) )
print(FeaturePlot(GW18GENormDat[[i]], features = c("ITM2A","IGFBP7"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW18GENormDat[[i]], features = c("MATN4","SCRG1"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW18GENormDat[[i]], features = c("SMAGP","HPGD"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW18GENormDat[[i]], features = c("ERMN","MOBP"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW18GENormDat[[i]], features = c("SLC44A1","GPR17"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW18GENormDat[[i]], features = c("C1QB","HEXB"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW18GENormDat[[i]], features = c("SLC17A6","GAD1"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW18GENormDat[[i]], features = c("AQP4","RAX"), pt.size = 0.2, ncol = 2))
dev.off()

## doublet id and removal -  
DefaultAssay(GW18GENormDat[[i]]) = 'SCT'
TMPsce = Seurat::as.SingleCellExperiment(GW18GENormDat[[i]])
TMPsce <- scDblFinder(TMPsce,clust.method='overcluster')
GW18GENormDat[[i]]@meta.data$scDblFinder = colData(TMPsce)$scDblFinder.class
GW18GENormDat[[i]] = subset(GW18GEAllDat[[i]],cells = colnames(GW18GENormDat[[i]])[which(GW18GENormDat[[i]]@meta.data$scDblFinder=='singlet')]) #7894
GW18GECellTracking[i,'AfterDoublets'] = ncol(GW18GENormDat[[i]])
#GW18GENormDat[[i]]@meta.data$Maj_Clust_From_Ref = AllRefGW18GE.integrated@meta.data[colnames(GW18GENormDat[[i]]),'Maj_Clust_From_Ref']
GW18GENormDat[[i]] <- SCTransform(GW18GENormDat[[i]], vars.to.regress = "percent.mt", verbose = FALSE)
GW18GENormDat[[i]] <- RunPCA(GW18GENormDat[[i]], verbose = FALSE)
GW18GENormDat[[i]] <- RunUMAP(GW18GENormDat[[i]], dims = 1:50, verbose = FALSE)
GW18GENormDat[[i]] <- FindNeighbors(GW18GENormDat[[i]], dims = 1:50, verbose = FALSE)
GW18GENormDat[[i]] <- FindClusters(GW18GENormDat[[i]], verbose = FALSE, resolution=0.5)

pdf(file=paste("./TestPlots/",i,"_Check_Cluster_PostDupRm.pdf",sep=''),width=12,height=8)
print(DimPlot(GW18GENormDat[[i]], label = TRUE) + NoLegend())
print(FeaturePlot(GW18GENormDat[[i]], features = c("AQP4","AGT"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW18GENormDat[[i]], features = c("FREM2","FZD2"), pt.size = 0.2, ncol = 2) )
print(FeaturePlot(GW18GENormDat[[i]], features = c("ITM2A","IGFBP7"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW18GENormDat[[i]], features = c("MATN4","SCRG1"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW18GENormDat[[i]], features = c("SMAGP","HPGD"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW18GENormDat[[i]], features = c("ERMN","MOBP"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW18GENormDat[[i]], features = c("SLC44A1","GPR17"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW18GENormDat[[i]], features = c("C1QB","HEXB"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW18GENormDat[[i]], features = c("SLC17A6","GAD1"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW18GENormDat[[i]], features = c("AQP4","RAX"), pt.size = 0.2, ncol = 2))
dev.off()

cat(paste('\n\n',i,' done normalization','\n\n',sep=''))

}

save(GW18GEAllDat,file='./SeuratObj/GW18GESamples_AllCells.rda',compress=TRUE)

## used in paper:
## load(url('https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/GW18GESamples_AllCells.rda'))

save(GW18GENormDat,file='./SeuratObj/GW18GESamples_PostSCT.rda',compress=TRUE)

## used in paper:
## load(url('https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/GW18GESamples_PostSCT.rda'))

write.csv(GW18GECellTracking,file='./SeuratObj/GW18GESamples_CellCounts.csv')


## Sample GW19 

GW19_Samples = c("GW19_CGE","GW19_LGE","GW19_MGE","GW19_motor","GW19_parietal","GW19_somato","GW19_PFC","GW19_V1")

GW19_motor_data <- Read10X(data.dir = "./Counts/GW19_M1_all" ) 
GW19_parietal_data <- Read10X(data.dir = "./Counts/GW19_Parietal")  
GW19_somato_data <- Read10X(data.dir = "./Counts/GW19_S1")
GW19_CGE_data <- Read10X(data.dir = "./Counts/GW19_CGE" )
GW19_LGE_data <- Read10X(data.dir = "./Counts/GW19_LGE"   )       
GW19_MGE_data <- Read10X(data.dir =  "./Counts/GW19_MGE" ) 
GW19_PFC_data <- Read10X(data.dir =  "./Counts/GW19_PFC_all")      
GW19_V1_data <- Read10X(data.dir =   "./Counts/GW19_V1_all" )

GW19_RawDat = list(GW19_CGE_data,GW19_LGE_data,GW19_MGE_data,GW19_motor_data,GW19_parietal_data,GW19_somato_data,GW19_PFC_data,GW19_V1_data)

names(GW19_RawDat) = GW19_Samples

setwd('/local/projects-t3/idea/bherb/Hypothalamus/PubRes/Rewrite')


GW19_CellTracking = data.frame(sample=GW19_Samples,Raw=0,AfterCutoffs=0,AfterDoublets=0)
rownames(GW19_CellTracking) = GW19_Samples

GW19_Dat = vector(mode='list',length=length(GW19_Samples))
names(GW19_Dat) = GW19_Samples
GW19_NormDat = vector(mode='list',length=length(GW19_Samples))
names(GW19_NormDat) = GW19_Samples


for(i in GW19_Samples){

GW19_Dat[[i]] <- CreateSeuratObject(counts = GW19_RawDat[[i]])
GW19_CellTracking[i,'Raw'] = ncol(GW19_Dat[[i]])
GW19_Dat[[i]]@meta.data$sample=i
GW19_Dat[[i]] = RenameCells(GW19_Dat[[i]],add.cell.id=i)
GW19_Dat[[i]] <- PercentageFeatureSet(GW19_Dat[[i]], pattern = "^MT-", col.name = "percent.mt")
## QC plots
pdf(file=paste("./TestPlots/",i,"_all_QC.pdf",sep=''),width=12,height=8)
print(VlnPlot(GW19_Dat[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
)
print(FeatureScatter(GW19_Dat[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt"))
print(FeatureScatter(GW19_Dat[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))
dev.off()
GW19_NormDat[[i]] = subset(GW19_Dat[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA > 1000& nCount_RNA < 15000 & percent.mt < 20)
GW19_CellTracking[i,'AfterCutoffs'] = ncol(GW19_NormDat[[i]])
if(ncol(GW19_NormDat[[i]])<100) next
GW19_NormDat[[i]] <- SCTransform(GW19_NormDat[[i]], vars.to.regress = "percent.mt", verbose = FALSE)
GW19_NormDat[[i]] <- RunPCA(GW19_NormDat[[i]], verbose = FALSE)
GW19_NormDat[[i]] <- RunUMAP(GW19_NormDat[[i]], dims = 1:50, verbose = FALSE)
GW19_NormDat[[i]] <- FindNeighbors(GW19_NormDat[[i]], dims = 1:50, verbose = FALSE)
GW19_NormDat[[i]] <- FindClusters(GW19_NormDat[[i]], verbose = FALSE, resolution=0.5)

DefaultAssay(GW19_NormDat[[i]]) = 'RNA'

pdf(file=paste("./TestPlots/",i,"_Check_Cluster_PreDupRm.pdf",sep=''),width=12,height=8)
print(DimPlot(GW19_NormDat[[i]], reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE))
print(FeaturePlot(GW19_NormDat[[i]], features = c("AQP4","AGT"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW19_NormDat[[i]], features = c("FREM2","FZD2"), pt.size = 0.2, ncol = 2) )
print(FeaturePlot(GW19_NormDat[[i]], features = c("ITM2A","IGFBP7"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW19_NormDat[[i]], features = c("MATN4","SCRG1"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW19_NormDat[[i]], features = c("SMAGP","HPGD"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW19_NormDat[[i]], features = c("ERMN","MOBP"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW19_NormDat[[i]], features = c("SLC44A1","GPR17"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW19_NormDat[[i]], features = c("C1QB","HEXB"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW19_NormDat[[i]], features = c("SLC17A6","GAD1"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW19_NormDat[[i]], features = c("AQP4","RAX"), pt.size = 0.2, ncol = 2))
dev.off()

## doublet id and removal -  
DefaultAssay(GW19_NormDat[[i]]) = 'SCT'
Seurat_Object_Diet <- DietSeurat(GW19_NormDat[[i]], graphs = "pca")
TMPsce <- as.SingleCellExperiment(Seurat_Object_Diet)
#TMPsce = as.SingleCellExperiment(GW19_NormDat[[i]])
TMPsce <- scDblFinder(TMPsce)
GW19_NormDat[[i]]@meta.data$scDblFinder = colData(TMPsce)$scDblFinder.class
GW19_NormDat[[i]] = subset(GW19_Dat[[i]],cells = colnames(GW19_NormDat[[i]])[which(GW19_NormDat[[i]]@meta.data$scDblFinder=='singlet')]) #7894
GW19_CellTracking[i,'AfterDoublets'] = ncol(GW19_NormDat[[i]])
if(ncol(GW19_NormDat[[i]])<100) next
GW19_NormDat[[i]] <- SCTransform(GW19_NormDat[[i]], vars.to.regress = "percent.mt", verbose = FALSE)
GW19_NormDat[[i]] <- RunPCA(GW19_NormDat[[i]], verbose = FALSE)
GW19_NormDat[[i]] <- RunUMAP(GW19_NormDat[[i]], dims = 1:50, verbose = FALSE)
GW19_NormDat[[i]] <- FindNeighbors(GW19_NormDat[[i]], dims = 1:50, verbose = FALSE)
GW19_NormDat[[i]] <- FindClusters(GW19_NormDat[[i]], verbose = FALSE, resolution=0.5)

pdf(file=paste("./TestPlots/",i,"_Check_Cluster_PostDupRm.pdf",sep=''),width=12,height=8)
print(DimPlot(GW19_NormDat[[i]], reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE))
print(FeaturePlot(GW19_NormDat[[i]], features = c("AQP4","AGT"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW19_NormDat[[i]], features = c("FREM2","FZD2"), pt.size = 0.2, ncol = 2) )
print(FeaturePlot(GW19_NormDat[[i]], features = c("ITM2A","IGFBP7"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW19_NormDat[[i]], features = c("MATN4","SCRG1"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW19_NormDat[[i]], features = c("SMAGP","HPGD"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW19_NormDat[[i]], features = c("ERMN","MOBP"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW19_NormDat[[i]], features = c("SLC44A1","GPR17"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW19_NormDat[[i]], features = c("C1QB","HEXB"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW19_NormDat[[i]], features = c("SLC17A6","GAD1"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW19_NormDat[[i]], features = c("AQP4","RAX"), pt.size = 0.2, ncol = 2))
dev.off()

cat(paste('\n\n',i,' done normalization','\n\n',sep=''))

}


save(GW19_Dat,file='./SeuratObj/GW19_NonHySamples_AllCells.rda',compress=TRUE)

## used in paper:
## load(url('https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/GW19_NonHySamples_AllCells.rda'))

save(GW19_NormDat,file='./SeuratObj/GW19_NonHySamples_PostSCT.rda',compress=TRUE)

## used in paper:
## load(url('https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/GW19_NonHySamples_AllCells.rda'))

write.csv(GW19_CellTracking,file='./SeuratObj/GW19_NonHySamples_CellCounts.csv')


## Sample GW20_34

GW20_34_Samples = c("GW20_34_CGE","GW20_34_LGE","GW20_34_MGE","GW20_34_parietal","GW20_34_PFC","GW20_34_V1")

GW20_34_parietal_data <- Read10X(data.dir = "./Counts/GW20_34_parietal")  
GW20_34_CGE_data <- Read10X(data.dir = "./Counts/GW20_34_CGE" )
GW20_34_LGE_data <- Read10X(data.dir = "./Counts/GW20_34_LGE"   )       
GW20_34_MGE_data <- Read10X(data.dir =  "./Counts/GW20_34_MGE" ) 
GW20_34_PFC_data <- Read10X(data.dir =  "./Counts/GW20_34_PFC"   )      
GW20_34_V1_data <- Read10X(data.dir =   "./Counts/GW20_34_V1" )

GW20_34_RawDat = list(GW20_34_CGE_data,GW20_34_LGE_data,GW20_34_MGE_data,GW20_34_parietal_data,GW20_34_PFC_data,GW20_34_V1_data)

names(GW20_34_RawDat) = GW20_34_Samples

GW20_34_CellTracking = data.frame(sample=GW20_34_Samples,Raw=0,AfterCutoffs=0,AfterDoublets=0)
rownames(GW20_34_CellTracking) = GW20_34_Samples

GW20_34_Dat = vector(mode='list',length=length(GW20_34_Samples))
names(GW20_34_Dat) = GW20_34_Samples
GW20_34_NormDat = vector(mode='list',length=length(GW20_34_Samples))
names(GW20_34_NormDat) = GW20_34_Samples


for(i in GW20_34_Samples){

GW20_34_Dat[[i]] <- CreateSeuratObject(counts = GW20_34_RawDat[[i]])
GW20_34_CellTracking[i,'Raw'] = ncol(GW20_34_Dat[[i]])
GW20_34_Dat[[i]]@meta.data$sample=i
GW20_34_Dat[[i]] = RenameCells(GW20_34_Dat[[i]],add.cell.id=i)
GW20_34_Dat[[i]] <- PercentageFeatureSet(GW20_34_Dat[[i]], pattern = "^MT-", col.name = "percent.mt")
## QC plots
pdf(file=paste("./TestPlots/",i,"_all_QC.pdf",sep=''),width=12,height=8)
print(VlnPlot(GW20_34_Dat[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
)
print(FeatureScatter(GW20_34_Dat[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt"))
print(FeatureScatter(GW20_34_Dat[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))
dev.off()
GW20_34_NormDat[[i]] = subset(GW20_34_Dat[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA > 1000& nCount_RNA < 15000 & percent.mt < 20)
GW20_34_CellTracking[i,'AfterCutoffs'] = ncol(GW20_34_NormDat[[i]])
if(ncol(GW20_34_NormDat[[i]])<100) next
GW20_34_NormDat[[i]] <- SCTransform(GW20_34_NormDat[[i]], vars.to.regress = "percent.mt", verbose = FALSE)
GW20_34_NormDat[[i]] <- RunPCA(GW20_34_NormDat[[i]], verbose = FALSE)
GW20_34_NormDat[[i]] <- RunUMAP(GW20_34_NormDat[[i]], dims = 1:50, verbose = FALSE)
GW20_34_NormDat[[i]] <- FindNeighbors(GW20_34_NormDat[[i]], dims = 1:50, verbose = FALSE)
GW20_34_NormDat[[i]] <- FindClusters(GW20_34_NormDat[[i]], verbose = FALSE, resolution=0.5)

DefaultAssay(GW20_34_NormDat[[i]]) = 'RNA'

pdf(file=paste("./TestPlots/",i,"_Check_Cluster_PreDupRm.pdf",sep=''),width=12,height=8)
print(DimPlot(GW20_34_NormDat[[i]], reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE))
print(FeaturePlot(GW20_34_NormDat[[i]], features = c("AQP4","AGT"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW20_34_NormDat[[i]], features = c("FREM2","FZD2"), pt.size = 0.2, ncol = 2) )
print(FeaturePlot(GW20_34_NormDat[[i]], features = c("ITM2A","IGFBP7"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW20_34_NormDat[[i]], features = c("MATN4","SCRG1"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW20_34_NormDat[[i]], features = c("SMAGP","HPGD"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW20_34_NormDat[[i]], features = c("ERMN","MOBP"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW20_34_NormDat[[i]], features = c("SLC44A1","GPR17"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW20_34_NormDat[[i]], features = c("C1QB","HEXB"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW20_34_NormDat[[i]], features = c("SLC17A6","GAD1"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW20_34_NormDat[[i]], features = c("AQP4","RAX"), pt.size = 0.2, ncol = 2))
dev.off()

## doublet id and removal -  
DefaultAssay(GW20_34_NormDat[[i]]) = 'SCT'
Seurat_Object_Diet <- DietSeurat(GW20_34_NormDat[[i]], graphs = "pca")
TMPsce <- as.SingleCellExperiment(Seurat_Object_Diet)
#TMPsce = as.SingleCellExperiment(GW20_34_NormDat[[i]])
TMPsce <- scDblFinder(TMPsce)
GW20_34_NormDat[[i]]@meta.data$scDblFinder = colData(TMPsce)$scDblFinder.class
GW20_34_NormDat[[i]] = subset(GW20_34_Dat[[i]],cells = colnames(GW20_34_NormDat[[i]])[which(GW20_34_NormDat[[i]]@meta.data$scDblFinder=='singlet')]) #7894
GW20_34_CellTracking[i,'AfterDoublets'] = ncol(GW20_34_NormDat[[i]])
if(ncol(GW20_34_NormDat[[i]])<100) next
GW20_34_NormDat[[i]] <- SCTransform(GW20_34_NormDat[[i]], vars.to.regress = "percent.mt", verbose = FALSE)
GW20_34_NormDat[[i]] <- RunPCA(GW20_34_NormDat[[i]], verbose = FALSE)
GW20_34_NormDat[[i]] <- RunUMAP(GW20_34_NormDat[[i]], dims = 1:50, verbose = FALSE)
GW20_34_NormDat[[i]] <- FindNeighbors(GW20_34_NormDat[[i]], dims = 1:50, verbose = FALSE)
GW20_34_NormDat[[i]] <- FindClusters(GW20_34_NormDat[[i]], verbose = FALSE, resolution=0.5)

pdf(file=paste("./TestPlots/",i,"_Check_Cluster_PostDupRm.pdf",sep=''),width=12,height=8)
print(DimPlot(GW20_34_NormDat[[i]], reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE))
print(FeaturePlot(GW20_34_NormDat[[i]], features = c("AQP4","AGT"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW20_34_NormDat[[i]], features = c("FREM2","FZD2"), pt.size = 0.2, ncol = 2) )
print(FeaturePlot(GW20_34_NormDat[[i]], features = c("ITM2A","IGFBP7"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW20_34_NormDat[[i]], features = c("MATN4","SCRG1"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW20_34_NormDat[[i]], features = c("SMAGP","HPGD"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW20_34_NormDat[[i]], features = c("ERMN","MOBP"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW20_34_NormDat[[i]], features = c("SLC44A1","GPR17"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW20_34_NormDat[[i]], features = c("C1QB","HEXB"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW20_34_NormDat[[i]], features = c("SLC17A6","GAD1"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(GW20_34_NormDat[[i]], features = c("AQP4","RAX"), pt.size = 0.2, ncol = 2))
dev.off()

cat(paste('\n\n',i,' done normalization','\n\n',sep=''))

}


save(GW20_34_Dat,file='./SeuratObj/GW20_34_NonHySamples_AllCells.rda',compress=TRUE)

## used in paper:
## load(url('https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/GW20_34_NonHySamples_AllCells.rda'))

save(GW20_34_NormDat,file='./SeuratObj/GW20_34_NonHySamples_PostSCT.rda',compress=TRUE)

## used in paper:
## load(url('https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/GW20_34_NonHySamples_PostSCT.rda'))

write.csv(GW20_34_CellTracking,file='./SeuratObj/GW20_34_NonHySamples_CellCounts.csv')



## Human Adult samples from Lein lab 


# 10X190-5  HTHtub = tuberal region of HTH
# 10X190-6 HTHtub = tuberal region of HTH
# 10X192-7 HTHso = supraoptic region of HTH
# 10X192-8 HTHso = supraoptic region of HTH
# 10X193-3 MN = mammillary nucleus
# 10X193-4 MN = mammillary nucleus
# 10X203-7 HTHpo = preoptic region of HTH
# 10X203-8 HTHpo = preoptic region of HTH


EdHypoSampleID = c("10X190-5","10X190-6","10X192-7","10X192-8","10X193-3","10X193-4","10X203-7","10X203-8")

EdHypoSamples = c("HTHtub_1","HTHtub_2","HTHso_1","HTHso_2","HTHmn_1","HTHmn_2","HTHpo_1","HTHpo_2")

HTHtub_1_data <- Read10X(data.dir = paste('./Counts/',EdHypoSampleID[1],sep=''))
HTHtub_2_data <- Read10X(data.dir = paste('./Counts/',EdHypoSampleID[2],sep=''))
HTHso_1_data <- Read10X(data.dir = paste('./Counts/',EdHypoSampleID[3],sep=''))
HTHso_2_data <- Read10X(data.dir = paste('./Counts/',EdHypoSampleID[4],sep=''))
HTHmn_1_data <- Read10X(data.dir = paste('./Counts/',EdHypoSampleID[5],sep=''))
HTHmn_2_data <- Read10X(data.dir = paste('./Counts/',EdHypoSampleID[6],sep=''))
HTHpo_1_data <- Read10X(data.dir = paste('./Counts/',EdHypoSampleID[7],sep=''))
HTHpo_2_data <- Read10X(data.dir = paste('./Counts/',EdHypoSampleID[8],sep=''))

EdHypoRawDat = list(HTHtub_1_data,HTHtub_2_data,HTHso_1_data,HTHso_2_data,HTHmn_1_data,HTHmn_2_data,HTHpo_1_data,HTHpo_2_data)
names(EdHypoRawDat) = EdHypoSamples

## cell number tracking 

EdHypoCellTracking = data.frame(sample=EdHypoSamples,Raw=0,AfterCutoffs=0,AfterDoublets=0)
rownames(EdHypoCellTracking) = EdHypoSamples

EdHypoAllDat = vector(mode='list',length=length(EdHypoSamples))
names(EdHypoAllDat) = EdHypoSamples
EdHypoNormDat = vector(mode='list',length=length(EdHypoSamples))
names(EdHypoNormDat) = EdHypoSamples

for(i in EdHypoSamples){

EdHypoAllDat[[i]] <- CreateSeuratObject(counts = EdHypoRawDat[[i]])
EdHypoCellTracking[i,'Raw'] = ncol(EdHypoAllDat[[i]])
EdHypoAllDat[[i]]@meta.data$sample=i
EdHypoAllDat[[i]] = RenameCells(EdHypoAllDat[[i]],add.cell.id=i)
EdHypoAllDat[[i]] <- PercentageFeatureSet(EdHypoAllDat[[i]], pattern = "^MT-", col.name = "percent.mt")
## QC plots
pdf(file=paste("./TestPlots/",i,"_all_QC.pdf",sep=''),width=12,height=8)
print(VlnPlot(EdHypoAllDat[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
)
print(FeatureScatter(EdHypoAllDat[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt"))
print(FeatureScatter(EdHypoAllDat[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))
dev.off()
EdHypoNormDat[[i]] = subset(EdHypoAllDat[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & nCount_RNA > 1000& nCount_RNA < 100000 & percent.mt < 20)
EdHypoCellTracking[i,'AfterCutoffs'] = ncol(EdHypoNormDat[[i]])
EdHypoNormDat[[i]] <- SCTransform(EdHypoNormDat[[i]])
EdHypoNormDat[[i]] <- RunPCA(EdHypoNormDat[[i]], verbose = FALSE)
EdHypoNormDat[[i]] <- RunUMAP(EdHypoNormDat[[i]], dims = 1:50, verbose = FALSE)
EdHypoNormDat[[i]] <- FindNeighbors(EdHypoNormDat[[i]], dims = 1:50, verbose = FALSE)
EdHypoNormDat[[i]] <- FindClusters(EdHypoNormDat[[i]], verbose = FALSE, resolution=0.5)

DefaultAssay(EdHypoNormDat[[i]]) = 'RNA'

pdf(file=paste("./TestPlots/",i,"_Check_Cluster_PreDupRm.pdf",sep=''),width=12,height=8)
print(DimPlot(EdHypoNormDat[[i]], reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE))
print(FeaturePlot(EdHypoNormDat[[i]], features = c("VIM","ASCL1"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(EdHypoNormDat[[i]], features = c("SNAP25","SLC32A1"), pt.size = 0.2, ncol = 2) )
print(FeaturePlot(EdHypoNormDat[[i]], features = c("SLC17A6","SLC17A8"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(EdHypoNormDat[[i]], features = c("LHX1","LHX5"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(EdHypoNormDat[[i]], features = c("LHX6","LHX9"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(EdHypoNormDat[[i]], features = c("GAD1","DLX1"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(EdHypoNormDat[[i]], features = c("DLX2","DLX5"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(EdHypoNormDat[[i]], features = c("ISL1","NR5A1"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(EdHypoNormDat[[i]], features = c("OTP","POU3F2"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(EdHypoNormDat[[i]], features = c("AQP4","RAX"), pt.size = 0.2, ncol = 2))
dev.off()

## doublet id and removal -  
DefaultAssay(EdHypoNormDat[[i]]) = 'SCT'
TMPsce = Seurat::as.SingleCellExperiment(EdHypoNormDat[[i]])
TMPsce <- scDblFinder(TMPsce,clust.method='overcluster')
EdHypoNormDat[[i]]@meta.data$scDblFinder = colData(TMPsce)$scDblFinder.class
EdHypoNormDat[[i]] = subset(EdHypoAllDat[[i]],cells = colnames(EdHypoNormDat[[i]])[which(EdHypoNormDat[[i]]@meta.data$scDblFinder=='singlet')]) #7894
EdHypoCellTracking[i,'AfterDoublets'] = ncol(EdHypoNormDat[[i]])
EdHypoNormDat[[i]] <- SCTransform(EdHypoNormDat[[i]])
EdHypoNormDat[[i]] <- RunPCA(EdHypoNormDat[[i]], verbose = FALSE)
EdHypoNormDat[[i]] <- RunUMAP(EdHypoNormDat[[i]], dims = 1:50, verbose = FALSE)
EdHypoNormDat[[i]] <- FindNeighbors(EdHypoNormDat[[i]], dims = 1:50, verbose = FALSE)
EdHypoNormDat[[i]] <- FindClusters(EdHypoNormDat[[i]], verbose = FALSE, resolution=0.5)

pdf(file=paste("./TestPlots/",i,"_Check_Cluster_PostDupRm.pdf",sep=''),width=12,height=8)
print(DimPlot(EdHypoNormDat[[i]], reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE))
print(FeaturePlot(EdHypoNormDat[[i]], features = c("VIM","ASCL1"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(EdHypoNormDat[[i]], features = c("SNAP25","SLC32A1"), pt.size = 0.2, ncol = 2) )
print(FeaturePlot(EdHypoNormDat[[i]], features = c("SLC17A6","SLC17A8"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(EdHypoNormDat[[i]], features = c("LHX1","LHX5"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(EdHypoNormDat[[i]], features = c("LHX6","LHX9"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(EdHypoNormDat[[i]], features = c("GAD1","DLX1"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(EdHypoNormDat[[i]], features = c("DLX2","DLX5"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(EdHypoNormDat[[i]], features = c("ISL1","NR5A1"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(EdHypoNormDat[[i]], features = c("OTP","POU3F2"), pt.size = 0.2, ncol = 2))
print(FeaturePlot(EdHypoNormDat[[i]], features = c("AQP4","RAX"), pt.size = 0.2, ncol = 2))
dev.off()

cat(paste('\n\n',i,' done normalization','\n\n',sep=''))

}

save(EdHypoAllDat,file='./SeuratObj/EdHypoSamples_AllCells.rda',compress=TRUE)

## used in paper:
## load(url('https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/EdHypoSamples_AllCells.rda'))

save(EdHypoNormDat,file='./SeuratObj/EdHypoSamples_PostSCT.rda',compress=TRUE)

## used in paper:
## load(url('https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/EdHypoSamples_PostSCT.rda'))

write.csv(EdHypoCellTracking,file='./SeuratObj/EdHypoSamples_CellCounts.csv')






