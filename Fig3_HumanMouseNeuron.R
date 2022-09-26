#R Studio Version 2022.07.1 Build 554
#base r v4.1.2
library(plotly) #v4.10.0
library(Seurat) #v4.1.1
library(tidyverse) #v1.3.2
library(plyr) #v1.8.7
library(htmlwidgets) #v1.5.4
library(monocle3) #v1.0.0
library(scales) #v1.2.1
library(igraph) #v1.3.4
library(ggraph) #v2.0.5
library(graphlayouts) #v0.8.0
library(qgraph) #v1.9.2
library(ggpubr) #v0.4.0

#### FIGURE 3: POMC CLustering

GeneLists = list("Arcuate" = c("SLC32A1","GAD1", "SLC17A6", "HDC", "POMC", "LEPR", "SST", "ISL1",  "SIX6", "SIX3", "HMX2", "GSX1", "KISS1", "GHRH", "TH", "TRH","BDNF", "PCSK1", "OTP", "ADCYAP1", "GAL","NPY","AGRP", "ASCL1", "HES1", "MKI67", "FABP7", "TTYH1", "HMGA2", "NKX2-1", "NHLH2","HOPX", "HES1", "NES", "SOX2", "RBFOX3", "NPY1R", "CRABP1"))


###### CLUSTERING POMC NEURONS
#load("~/Downloads/AllHSMMNeuronsPure_integrated.rda")
load(url('https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/AllHSMMNeuronsPure_integrated.rda')) # All.integrated
DefaultAssay(All.integrated) = "RNA"
FeaturePlot(All.integrated, "POMC")

Idents(All.integrated) ="seurat_clusters"
DimPlot(All.integrated, label=T)
All.integrated_POMC = subset(All.integrated, idents = c(20))


DefaultAssay(All.integrated_POMC) = "integrated"
All.integrated_POMC = RunPCA(All.integrated_POMC, npcs = 10)
All.integrated_POMC <- RunUMAP(All.integrated_POMC, dims = 1:10, spread= 2)


set.dim = 10
set.kparam = 10
set.res = 4
ClusterFunc_All_RNA(All.integrated_POMC)

DefaultAssay(All.integrated_POMC) = "integrated"
All.integrated_POMC <- FindNeighbors(All.integrated_POMC, k.param=50, dims=1:10)
All.integrated_POMC <- FindClusters(All.integrated_POMC, resolution = 2)
All.integrated_POMC_Clean = subset(All.integrated_POMC, idents = c(11), invert=T)


#### ITERATE CLUSTERS
DPlist2 = list()
DPlist = list()
FPlist = list()
Filename = "All.integrated_POMC_Clean_17AUG22"
for(y in seq(10,50,10)){
  DefaultAssay(All.integrated_POMC_Clean) = "integrated"
  All.integrated_POMC_Clean = RunPCA(All.integrated_POMC_Clean, npcs = y)
  for(z in seq(10,y,10)){
    for(s in c(2,5,10)){
      All.integrated_POMC_Clean <- RunUMAP(All.integrated_POMC_Clean, dims = 1:z, spread= s)
      DefaultAssay(All.integrated_POMC_Clean) = "RNA"
      
      FPlist[[paste("PCA", y, "_dims", z, "_spread", s)]] = FeaturePlot(All.integrated_POMC_Clean, c("SLC32A1","GAD1", "SLC17A6", "HDC", "POMC", "LEPR", "SST", "ISL1",  "SIX6", "SIX3", "HMX2", "GSX1", "KISS1", "GHRH", "TH", "TRH","BDNF", "PCSK1", "OTP", "ADCYAP1", "GAL","NPY","AGRP", "ASCL1", "HES1", "MKI67", "FABP7", "TTYH1", "HMGA2", "NKX2-1", "NHLH2","HOPX", "HES1", "NES", "SOX2", "RBFOX3", "NPY1R", "CRABP1"), reduction="umap") 
      
      DPlist[[paste("PCA", y, "_dims", z, "_spread", s)]] = DimPlot(All.integrated_POMC_Clean, reduction="umap", split.by = "sample") + labs(title = paste("PCA", y, "_dims", z, "_spread", s))
      
      DPlist2[[paste("PCA", y, "_dims", z, "_spread", s)]] = DimPlot(All.integrated_POMC_Clean, reduction="umap", group.by = "sample", label=T) + labs(title = paste("PCA", y, "dims", z, "spread", s))
    }}}  
pdf(paste("", Filename, "_FeaturePlots.pdf", sep=""), width=20, height=40)
print(FPlist)
dev.off()

pdf(paste("", Filename, "_SPLITUMAP.pdf", sep=""), width=50, height=5)
print(DPlist)
dev.off()

pdf(paste("", Filename, "_GROUPUMAP.pdf", sep=""), width=8, height=5)
print(DPlist2)
dev.off()
####

DefaultAssay(All.integrated_POMC_Clean) = "integrated"
All.integrated_POMC_Clean = RunPCA(All.integrated_POMC_Clean, npcs = 10)
All.integrated_POMC_Clean <- RunUMAP(All.integrated_POMC_Clean, dims = 1:10, spread= 2)


set.dim = 10
set.kparam = 10
set.res = 1
ClusterFunc_All_RNA(All.integrated_POMC_Clean)

DefaultAssay(All.integrated_POMC_Clean) = "integrated"
All.integrated_POMC_Clean <- FindNeighbors(All.integrated_POMC_Clean, k.param=10, dims=1:10)
All.integrated_POMC_Clean <- FindClusters(All.integrated_POMC_Clean, resolution = 1)
All.integrated_POMC_Clean_UMAP = as.data.frame(All.integrated_POMC_Clean@reductions$umap@cell.embeddings)
POMC_01 = subset(All.integrated_POMC_Clean, idents = c(0)) #POMC, ISL1, SIX3, NPY1R
POMC_02 = subset(All.integrated_POMC_Clean, idents = c(1,3,11)) #POMC, GAD1, SLC17A6, ISL1, SIX3, PCSK1, NPY1R
POMC_03 = subset(All.integrated_POMC_Clean, idents = c(2)) #POMC, GAD1, SLC17A6, ISL1, SIX6, SIX3, PCSK1, NKX2-1, NPY1R
POMC_04 = subset(All.integrated_POMC_Clean, idents = c(4)) #POMC, ISL1, SIX6, SIX3, NKX2-1, 
POMC_05 = subset(All.integrated_POMC_Clean, idents = c(6, 13)) #POMC, ISL1, SIX6, SIX3, NKX2-1, 
POMC_06 = subset(All.integrated_POMC_Clean, idents = c(7)) #POMC, GAD1, SLC17A6, LEPR, ISL1, SIX6, SIX3, PCSK1, NKX2-1, NPY1R
POMC_07 = subset(All.integrated_POMC_Clean, idents = c(8, 10)) #POMC, GAD1, SLC17A6, ISL1, SIX6, SIX3, PCSK1, NPY1R
POMC_08 = subset(All.integrated_POMC_Clean, idents = c(9)) #POMC, ISL1, SIX3, NKX2-1, 
POMC_09 = subset(All.integrated_POMC_Clean, idents = c(5, 12, 16)) #POMC, ISL1, SIX6, SIX3, NKX2-1
POMC_10 = subset(All.integrated_POMC_Clean, idents = c(14)) #POMC, SLC17A6, SST, ISL1, SIX6, SIX3, PCSK1, NPY1R
POMC_11 = subset(All.integrated_POMC_Clean, idents = c(15)) #POMC, SLC17A6, LEPR, ISL1, SIX3, PCSK1, NPY1R

#CLEANING STRAYS
POMC_02_to06 = subset(All.integrated_POMC_Clean_UMAP, row.names(All.integrated_POMC_Clean_UMAP) %in% colnames(POMC_02) &  All.integrated_POMC_Clean_UMAP$UMAP_1 > -4 & All.integrated_POMC_Clean_UMAP$UMAP_1 < 0 &  All.integrated_POMC_Clean_UMAP$UMAP_2 > 2 & All.integrated_POMC_Clean_UMAP$UMAP_2 < 5 )
POMC_02_to07 = subset(All.integrated_POMC_Clean_UMAP, row.names(All.integrated_POMC_Clean_UMAP) %in% colnames(POMC_02) &  All.integrated_POMC_Clean_UMAP$UMAP_1 > -4 & All.integrated_POMC_Clean_UMAP$UMAP_1 < 0 &  All.integrated_POMC_Clean_UMAP$UMAP_2 < 2 | row.names(All.integrated_POMC_Clean_UMAP) %in% colnames(POMC_02) &  All.integrated_POMC_Clean_UMAP$UMAP_1 > -7 & All.integrated_POMC_Clean_UMAP$UMAP_1 < 0 &  All.integrated_POMC_Clean_UMAP$UMAP_2 < -1.5)
POMC_02_to01 = subset(All.integrated_POMC_Clean_UMAP, row.names(All.integrated_POMC_Clean_UMAP) %in% colnames(POMC_02) &  All.integrated_POMC_Clean_UMAP$UMAP_1 >  0 &  All.integrated_POMC_Clean_UMAP$UMAP_2 > 0)
POMC_02_to03 = subset(All.integrated_POMC_Clean_UMAP, row.names(All.integrated_POMC_Clean_UMAP) %in% colnames(POMC_02) &  All.integrated_POMC_Clean_UMAP$UMAP_1 >  0 &  All.integrated_POMC_Clean_UMAP$UMAP_2 < 0)
POMC_2_Rest = subset(POMC_02, cells = c(row.names(POMC_02_to06), row.names(POMC_02_to07), row.names(POMC_02_to01), row.names(POMC_02_to03)), invert=T)

POMC_03_to05 = subset(All.integrated_POMC_Clean_UMAP, row.names(All.integrated_POMC_Clean_UMAP) %in% colnames(POMC_03) &  All.integrated_POMC_Clean_UMAP$UMAP_2 > 5)
POMC_3_Rest = subset(POMC_03, cells = c(row.names(POMC_03_to05)), invert=T)

POMC_04_to06 = subset(All.integrated_POMC_Clean_UMAP, row.names(All.integrated_POMC_Clean_UMAP) %in% colnames(POMC_04) &  All.integrated_POMC_Clean_UMAP$UMAP_1 > -4 & All.integrated_POMC_Clean_UMAP$UMAP_1 < 0 &  All.integrated_POMC_Clean_UMAP$UMAP_2 > 2 & All.integrated_POMC_Clean_UMAP$UMAP_2 < 5 )
POMC_04_to07 = subset(All.integrated_POMC_Clean_UMAP, row.names(All.integrated_POMC_Clean_UMAP) %in% colnames(POMC_04) &  All.integrated_POMC_Clean_UMAP$UMAP_1 > -4 & All.integrated_POMC_Clean_UMAP$UMAP_1 < 0 &  All.integrated_POMC_Clean_UMAP$UMAP_2 < 2)
POMC_04_to09 = subset(All.integrated_POMC_Clean_UMAP, row.names(All.integrated_POMC_Clean_UMAP) %in% colnames(POMC_04) &  All.integrated_POMC_Clean_UMAP$UMAP_1 > 2 &  All.integrated_POMC_Clean_UMAP$UMAP_2 > 1)
POMC_04_to01 = subset(All.integrated_POMC_Clean_UMAP, row.names(All.integrated_POMC_Clean_UMAP) %in% colnames(POMC_04) &  All.integrated_POMC_Clean_UMAP$UMAP_1 > 0 & All.integrated_POMC_Clean_UMAP$UMAP_1 < 2 &  All.integrated_POMC_Clean_UMAP$UMAP_2 > 1)
POMC_4_Rest = subset(POMC_04, cells = c(row.names(POMC_04_to06), row.names(POMC_04_to07), row.names(POMC_04_to09), row.names(POMC_04_to01)), invert=T)

POMC_05_to10 = subset(All.integrated_POMC_Clean_UMAP, row.names(All.integrated_POMC_Clean_UMAP) %in% colnames(POMC_05) &   All.integrated_POMC_Clean_UMAP$UMAP_1 < 0 & All.integrated_POMC_Clean_UMAP$UMAP_2 > 4)
POMC_05_to01 = subset(All.integrated_POMC_Clean_UMAP, row.names(All.integrated_POMC_Clean_UMAP) %in% colnames(POMC_05) &  All.integrated_POMC_Clean_UMAP$UMAP_1 > -2 & All.integrated_POMC_Clean_UMAP$UMAP_2 < 2.5 )
POMC_5_Rest = subset(POMC_05, cells = c(row.names(POMC_05_to10), row.names(POMC_05_to01)), invert=T)

POMC_07_to01 = subset(All.integrated_POMC_Clean_UMAP, row.names(All.integrated_POMC_Clean_UMAP) %in% colnames(POMC_07) &  All.integrated_POMC_Clean_UMAP$UMAP_1 > -1 & All.integrated_POMC_Clean_UMAP$UMAP_1 < 3 & All.integrated_POMC_Clean_UMAP$UMAP_2 < 4)
POMC_07_to02 = subset(All.integrated_POMC_Clean_UMAP, row.names(All.integrated_POMC_Clean_UMAP) %in% colnames(POMC_07) &  All.integrated_POMC_Clean_UMAP$UMAP_1 < -6)
POMC_07_to06 = subset(All.integrated_POMC_Clean_UMAP, row.names(All.integrated_POMC_Clean_UMAP) %in% colnames(POMC_07) &  All.integrated_POMC_Clean_UMAP$UMAP_1 > -4 & All.integrated_POMC_Clean_UMAP$UMAP_1 < 0 &  All.integrated_POMC_Clean_UMAP$UMAP_2 > 2 & All.integrated_POMC_Clean_UMAP$UMAP_2 < 5 )
POMC_07_to05 = subset(All.integrated_POMC_Clean_UMAP, row.names(All.integrated_POMC_Clean_UMAP) %in% colnames(POMC_07) &  All.integrated_POMC_Clean_UMAP$UMAP_2 > 5)
POMC_07_to09 = subset(All.integrated_POMC_Clean_UMAP, row.names(All.integrated_POMC_Clean_UMAP) %in% colnames(POMC_07) &  All.integrated_POMC_Clean_UMAP$UMAP_1 > 4)
POMC_7_Rest = subset(POMC_07, cells = c(row.names(POMC_07_to02), row.names(POMC_07_to06), row.names(POMC_07_to05), row.names(POMC_07_to09),  row.names(POMC_07_to01)), invert=T)


POMC_08_to01 = subset(All.integrated_POMC_Clean_UMAP, row.names(All.integrated_POMC_Clean_UMAP) %in% colnames(POMC_08) &  All.integrated_POMC_Clean_UMAP$UMAP_1 > -1 & All.integrated_POMC_Clean_UMAP$UMAP_2 > -3)
POMC_08_to07 = subset(All.integrated_POMC_Clean_UMAP, row.names(All.integrated_POMC_Clean_UMAP) %in% colnames(POMC_08) & All.integrated_POMC_Clean_UMAP$UMAP_1 < 0)
POMC_08_to03 = subset(All.integrated_POMC_Clean_UMAP, row.names(All.integrated_POMC_Clean_UMAP) %in% colnames(POMC_08) &  All.integrated_POMC_Clean_UMAP$UMAP_1 < 4.5 & All.integrated_POMC_Clean_UMAP$UMAP_1 > 0 & All.integrated_POMC_Clean_UMAP$UMAP_2 < -3)
POMC_8_Rest = subset(POMC_08, cells = c(row.names(POMC_08_to01), row.names(POMC_08_to07), row.names(POMC_08_to03)), invert=T)

POMC_09_to07 = subset(All.integrated_POMC_Clean_UMAP, row.names(All.integrated_POMC_Clean_UMAP) %in% colnames(POMC_09) & All.integrated_POMC_Clean_UMAP$UMAP_1 < 0 & All.integrated_POMC_Clean_UMAP$UMAP_2 < 0)
POMC_09_to06 = subset(All.integrated_POMC_Clean_UMAP, row.names(All.integrated_POMC_Clean_UMAP) %in% colnames(POMC_09) &  All.integrated_POMC_Clean_UMAP$UMAP_1 > -4 & All.integrated_POMC_Clean_UMAP$UMAP_1 < 0 &  All.integrated_POMC_Clean_UMAP$UMAP_2 > 2 & All.integrated_POMC_Clean_UMAP$UMAP_2 < 5 )
POMC_09_to10 = subset(All.integrated_POMC_Clean_UMAP, row.names(All.integrated_POMC_Clean_UMAP) %in% colnames(POMC_09) &   All.integrated_POMC_Clean_UMAP$UMAP_1 < 0 & All.integrated_POMC_Clean_UMAP$UMAP_2 > 4)
POMC_09_to01 = subset(All.integrated_POMC_Clean_UMAP, row.names(All.integrated_POMC_Clean_UMAP) %in% colnames(POMC_09) &  All.integrated_POMC_Clean_UMAP$UMAP_1 > 0  & All.integrated_POMC_Clean_UMAP$UMAP_1 < 3.5 & All.integrated_POMC_Clean_UMAP$UMAP_2 < 4)
POMC_09_to05 = subset(All.integrated_POMC_Clean_UMAP, row.names(All.integrated_POMC_Clean_UMAP) %in% colnames(POMC_09) &  All.integrated_POMC_Clean_UMAP$UMAP_2 > 5 & All.integrated_POMC_Clean_UMAP$UMAP_1 > 2 & All.integrated_POMC_Clean_UMAP$UMAP_1 < 5)
POMC_9_Rest = subset(POMC_09, cells = c(row.names(POMC_09_to07), row.names(POMC_09_to06), row.names(POMC_09_to10), row.names(POMC_09_to01), row.names(POMC_09_to05)), invert=T)

POMC_10_to02 = subset(All.integrated_POMC_Clean_UMAP, row.names(All.integrated_POMC_Clean_UMAP) %in% colnames(POMC_10) &  All.integrated_POMC_Clean_UMAP$UMAP_1 < -5 )
POMC_10_to01 = subset(All.integrated_POMC_Clean_UMAP, row.names(All.integrated_POMC_Clean_UMAP) %in% colnames(POMC_10) &  All.integrated_POMC_Clean_UMAP$UMAP_1 > -5 & All.integrated_POMC_Clean_UMAP$UMAP_2 < 4)
POMC_10_Rest = subset(POMC_10, cells = c(row.names(POMC_10_to02), row.names(POMC_10_to01)), invert=T)


POMC_11_to02 = subset(All.integrated_POMC_Clean_UMAP, row.names(All.integrated_POMC_Clean_UMAP) %in% colnames(POMC_11) &  All.integrated_POMC_Clean_UMAP$UMAP_1 < -5 & All.integrated_POMC_Clean_UMAP$UMAP_2 < 2 )
POMC_11_to05 = subset(All.integrated_POMC_Clean_UMAP, row.names(All.integrated_POMC_Clean_UMAP) %in% colnames(POMC_11) &  All.integrated_POMC_Clean_UMAP$UMAP_2 > 5)
POMC_11_to01 = subset(All.integrated_POMC_Clean_UMAP, row.names(All.integrated_POMC_Clean_UMAP) %in% colnames(POMC_11) &  All.integrated_POMC_Clean_UMAP$UMAP_1 > 0 & All.integrated_POMC_Clean_UMAP$UMAP_2 < 4)
POMC_11_Rest = subset(POMC_11, cells = c(row.names(POMC_11_to02), row.names(POMC_11_to05), row.names(POMC_11_to01)), invert=T)


###CleanedVersions
POMC_01_Cleaned = subset(All.integrated_POMC_Clean, cells = c(colnames(POMC_01), row.names(POMC_02_to01), row.names(POMC_04_to01), row.names(POMC_05_to01), row.names(POMC_08_to01), row.names(POMC_09_to01), row.names(POMC_10_to01), row.names(POMC_11_to01), row.names(POMC_07_to01)))
POMC_02_Cleaned = subset(All.integrated_POMC_Clean, cells = c(colnames(POMC_2_Rest), row.names(POMC_07_to02), row.names(POMC_10_to02), row.names(POMC_11_to02)))
POMC_03_Cleaned = subset(All.integrated_POMC_Clean, cells = c(colnames(POMC_3_Rest), row.names(POMC_02_to03), row.names(POMC_08_to03)))
POMC_04_Cleaned = subset(All.integrated_POMC_Clean, cells = c(colnames(POMC_4_Rest)))
POMC_05_Cleaned = subset(All.integrated_POMC_Clean, cells = c(colnames(POMC_5_Rest), row.names(POMC_03_to05), row.names(POMC_07_to05), row.names(POMC_11_to05)))
POMC_06_Cleaned = subset(All.integrated_POMC_Clean, cells = c(colnames(POMC_06), row.names(POMC_02_to06), row.names(POMC_04_to06), row.names(POMC_07_to06), row.names(POMC_09_to06)))
POMC_07_Cleaned = subset(All.integrated_POMC_Clean, cells = c(colnames(POMC_7_Rest), row.names(POMC_02_to07), row.names(POMC_04_to07), row.names(POMC_08_to07), row.names(POMC_09_to07)))
POMC_08_Cleaned = subset(All.integrated_POMC_Clean, cells = c(colnames(POMC_8_Rest)))
POMC_09_Cleaned = subset(All.integrated_POMC_Clean, cells = c(colnames(POMC_9_Rest), row.names(POMC_04_to09), row.names(POMC_07_to09), row.names(POMC_09_to05)))
POMC_10_Cleaned = subset(All.integrated_POMC_Clean, cells = c(colnames(POMC_10_Rest), row.names(POMC_05_to10), row.names(POMC_09_to10)))
POMC_11_Cleaned = subset(All.integrated_POMC_Clean, cells = c(colnames(POMC_11_Rest)))


POMC_Barcs = GenerateMetaData(list("POMC_01" = POMC_01_Cleaned, "POMC_02" = POMC_02_Cleaned,  "POMC_04" = POMC_04_Cleaned, "POMC_05" = POMC_05_Cleaned, "POMC_06" = POMC_06_Cleaned, "POMC_07" = POMC_07_Cleaned, "POMC_08" = POMC_08_Cleaned, "POMC_09" = POMC_09_Cleaned, "POMC_10" = POMC_10_Cleaned, "POMC_11" = POMC_11_Cleaned, "POMC_03" = POMC_03_Cleaned))


#### FIGURE 3

load("~/Downloads/AllHSMMNeuronsPure_integrated.rda")
DefaultAssay(All.integrated) = "RNA"
All.integrated@meta.data$sampletype = gsub("HTHtub_1", "HumanAdult - Tub", gsub("HTHtub_2", "HumanAdult - Tub", gsub("HTHso_1", "HumanAdult - SO", gsub("HTHso_2", "HumanAdult - SO", gsub("HTHmn_1", "HumanAdult - Mammill", gsub("HTHmn_2", "HumanAdult - Mammill", gsub("HTHpo_1", "HumanAdult - PO", gsub("HTHpo_2", "HumanAdult - PO", gsub("CS22_hypo", "HumanFetal - CS22", gsub("CS22_2_hypo", "HumanFetal - CS22", gsub("GW16_hypo", "HumanFetal - GW16", gsub("GW18_hypo", "HumanFetal - GW18", gsub("GW19_hypo", "HumanFetal - GW19", gsub("GW20_34_hypo", "HumanFetal - GW20", gsub("GW22T_hypo1", "HumanFetal - GW22", gsub("GW25_3V_hypo", "HumanFetal - GW25", gsub("Mof", "MouseAdult - PO", gsub("Chen", "MouseAdult", gsub("Mick.male", "MouseAdult - LH", gsub("Mick.female", "MouseAdult - LH", gsub("Kim_FC1", "MouseAdult - VMH", gsub("Kim_MC1", "MouseAdult - VMH", gsub("Kim_MC2", "MouseAdult - VMH", gsub("Kim_MC3", "MouseAdult - VMH", gsub("Wen", "MouseAdult - SCN", gsub("Cam", "MouseAdult - ARC", gsub("E10", "MouseFetal - E10", gsub("E11", "MouseFetal - E12", gsub("E12", "MouseFetal - E12", gsub("E13", "MouseFetal - E13", gsub("E14", "MouseFetal - E14", gsub("E15L", "MouseFetal - E15", gsub("\\bE15\\b" , "MouseFetal - E15", gsub("E16", "MouseFetal - E16", gsub("E16V1", "MouseFetal - E16", gsub("E18", "MouseFetal - E18", gsub("P4", "MouseAdult - P4", gsub("P8", "MouseAdult - P8", gsub("P14", "MouseAdult - P14", gsub("P45X", "P45", All.integrated@meta.data$sample))))))))))))))))))))))))))))))))))))))))
All.integrated@meta.data$sampletype = factor(All.integrated@meta.data$sampletype, levels = c("HumanFetal - CS22", "HumanFetal - GW16", "HumanFetal - GW18", "HumanFetal - GW19", "HumanFetal - GW20", "HumanFetal - GW22", "HumanFetal - GW25", "HumanAdult - Mammill", "HumanAdult - PO", "HumanAdult - SO", "HumanAdult - Tub", "MouseFetal - E10", "MouseFetal - E12", "MouseFetal - E13", "MouseFetal - E14", "MouseFetal - E15", "MouseFetal - E16", "MouseFetal - E18", "MouseAdult",    "MouseAdult - ARC",   "MouseAdult - PO",   "MouseAdult - SCN",  "MouseAdult - VMH", "MouseAdult - LH",
"MouseAdult - P4","MouseAdult - P8", "MouseAdult - P14",  "MouseAdult - P45"))
All.integrated@meta.data$sampletypeBroad = gsub(" -.*", "", All.integrated@meta.data$sampletype)

titlesz = 10
btmmarg = 0
Fig3List = list()
Idents(All.integrated) = "sampletype"
############## FIGURE 3A ##############
Fig3List[["AllDataUMAP"]] = DimPlot(All.integrated, raster=F, pt.size = 0.1)+NoLegend()

############## FIGURE 3B ##############
Idents(All.integrated) = "sampletypeBroad"
Fig3List[["SplitTimeUMAP"]] = DimPlot(All.integrated, raster=F, split.by = "sampletypeBroad", ncol =2, pt.size = 0.01)+NoLegend()+ theme(plot.title = element_text(face = "bold", size = titlesz, margin = margin(0,0,btmmarg,0)))

############## FIGURE 3C ##############
Fig3List[["AllDataFP"]] = FeaturePlot(All.integrated, c("SLC32A1","GAD1", "SLC17A6", "HDC"), raster=F, ncol =2, pt.size = 0.01)+ theme(plot.title = element_text(face = "bold", size = titlesz, margin = margin(0,0,btmmarg,0)))
Fig3List[["AllDataFP"]][[1]] = Fig3List[["AllDataFP"]][[1]] + theme(plot.title = element_text(face = "bold", size = titlesz, margin = margin(0,0,btmmarg,0)))
Fig3List[["AllDataFP"]][[2]] = Fig3List[["AllDataFP"]][[2]] + theme(plot.title = element_text(face = "bold", size = titlesz, margin = margin(0,0,btmmarg,0)))
Fig3List[["AllDataFP"]][[3]] = Fig3List[["AllDataFP"]][[3]] + theme(plot.title = element_text(face = "bold", size = titlesz, margin = margin(0,0,btmmarg,0)))
Fig3List[["AllDataFP"]][[4]] = Fig3List[["AllDataFP"]][[4]] + theme(plot.title = element_text(face = "bold", size = titlesz, margin = margin(0,0,btmmarg,0)))

############## FIGURE 3D ##############
Idents(All.integrated) = "seurat_clusters"
Fig3List[["AllData_Seu"]] = DimPlot(All.integrated, raster=F, pt.size = 0.1, label=T)+NoLegend()

############## FIGURE 3E ##############
All.integrated_nsamples = as.data.frame(table(All.integrated@meta.data$sampletype))

CompilePercents = as.data.frame(matrix(ncol =7, nrow=0))
colnames(CompilePercents) = c("SampleType", "nCluster", "nTotal", "PercentAll", "PercentClust_All",  "Cluster")
for(x in unique(All.integrated@meta.data$seurat_clusters)){
Subs= subset(All.integrated@meta.data, All.integrated@meta.data$seurat_clusters == x) 
MetaTable = as.data.frame(table(Subs$sampletype))
MetaTable2 = merge(MetaTable, All.integrated_nsamples, by = "Var1", all = T)
MetaTable2[is.na(MetaTable2)] = 0
colnames(MetaTable2) = c("SampleType", "nCluster", "nTotal")
MetaTable2$nCluster = as.numeric(MetaTable2$nCluster)
MetaTable2$nTotal = as.numeric(MetaTable2$nTotal)
MetaTable2[is.na(MetaTable2)] = 0
MetaTable2$PercentAll = MetaTable2$nCluster/MetaTable2$nTotal
MetaTable2$PercentClust_All = MetaTable2$PercentAll/sum(MetaTable2$PercentAll)*100
MetaTable2[is.na(MetaTable2)] = 0
MetaTable2$Cluster = x
CompilePercents = rbind(CompilePercents, MetaTable2)
}

CompilePercents$SampleType = factor(CompilePercents$SampleType, levels = c("HumanFetal - CS22", "HumanFetal - GW16", "HumanFetal - GW18", "HumanFetal - GW19", "HumanFetal - GW20", "HumanFetal - GW22", "HumanFetal - GW25", "HumanAdult - Mammill", "HumanAdult - PO", "HumanAdult - SO", "HumanAdult - Tub", "MouseFetal - E10", "MouseFetal - E12", "MouseFetal - E13", "MouseFetal - E14", "MouseFetal - E15", "MouseFetal - E16", "MouseFetal - E18", "MouseAdult",    "MouseAdult - ARC",   "MouseAdult - PO",   "MouseAdult - SCN",  "MouseAdult - VMH", "MouseAdult - LH",
"MouseAdult - P4","MouseAdult - P8", "MouseAdult - P14",  "MouseAdult - P45"))
CompilePercents$Cluster = factor(CompilePercents$Cluster, levels = sort(unique(All.integrated@meta.data$seurat_clusters)))
  
Fig3List[["Seu_NormalizedPercent"]] = ggplot(CompilePercents, aes(fill=SampleType, y=PercentClust_All, x=Cluster)) + geom_bar(position="stack", stat="identity") + ylab("% cells") + xlab("") + scale_y_continuous( expand= c(0,0)) + theme_classic() + theme(legend.position = "none", axis.text = element_text(size = 16), axis.title.y = element_text(size = 16), legend.title = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5)) 

############## FIGURE 3F ##############
Fig3List[["POMC_FP"]] = FeaturePlot(All.integrated, c("POMC"), raster=F, pt.size = 0.01)+ theme(plot.title = element_text(face = "bold", size = titlesz, margin = margin(0,0,btmmarg,0)))

############## FIGURE 3G ##############
All.integrated_POMC_Clean@meta.data$POMC_BarcsReorder = gsub("POMC_01", "POMC_B", gsub("POMC_02", "POMC_A", gsub("POMC_03", "POMC_D", gsub("POMC_04", "POMC_F", gsub("POMC_06", "POMC_H", gsub("POMC_08", "POMC_I", gsub("POMC_09", "POMC_C", All.integrated_POMC_Clean@meta.data$POMC_Barcs)))))))
All.integrated_POMC_Clean@meta.data$POMC_BarcsReorder = gsub("POMC_B", "POMC_02", gsub("POMC_A", "POMC_01", gsub("POMC_D", "POMC_04", gsub("POMC_F", "POMC_06", gsub("POMC_H", "POMC_08", gsub("POMC_I", "POMC_09", gsub("POMC_C", "POMC_03", All.integrated_POMC_Clean@meta.data$POMC_BarcsReorder)))))))
All.integrated_POMC_Clean@meta.data$POMC_BarcsReorder = factor(All.integrated_POMC_Clean@meta.data$POMC_BarcsReorder, levels = sort(unique(All.integrated_POMC_Clean@meta.data$POMC_BarcsReorder)))

Idents(All.integrated_POMC_Clean) = "POMC_BarcsReorder"
Fig3List[["POMC_Clustering"]] = DimPlot(All.integrated_POMC_Clean)+NoLegend()

############## FIGURE 3H ##############
POMCAllMeta = All.integrated_POMC_Clean@meta.data

CompilePercents = as.data.frame(matrix(ncol =9, nrow=0))
colnames(CompilePercents) = c("SampleType", "nCluster", "nTotal", "nPOMCneurons", "PercentAll", "PercentClust_All", "PercentPOMC", "PercentClust_POMC", "Timepoint")
for(x in unique(POMCAllMeta$POMC_BarcsReorder)){
Subs= subset(POMCAllMeta, POMCAllMeta$POMC_BarcsReorder == x) 
MetaTable = as.data.frame(table(Subs$sampletype))
MetaTable2 = merge(MetaTable, All.integrated_nsamples, by = "Var1", all = T)
MetaTable3 = merge(MetaTable2, All.integrated_POMC_Clean_nsamples, by = "Var1", all = T)
MetaTable3[is.na(MetaTable3)] = 0
colnames(MetaTable3) = c("SampleType", "nCluster", "nTotal", "nPOMCneurons")
MetaTable3$nCluster = as.numeric(MetaTable3$nCluster)
MetaTable3$nTotal = as.numeric(MetaTable3$nTotal)
MetaTable3$nPOMCneurons = as.numeric(MetaTable3$nPOMCneurons)
MetaTable3[is.na(MetaTable3)] = 0
MetaTable3$PercentAll = MetaTable3$nCluster/MetaTable3$nTotal
MetaTable3$PercentClust_All = MetaTable3$PercentAll/sum(MetaTable3$PercentAll)*100
MetaTable3$PercentPOMC = MetaTable3$nCluster/MetaTable3$nPOMCneurons
MetaTable3[is.na(MetaTable3)] = 0
MetaTable3$PercentClust_POMC = MetaTable3$PercentPOMC/sum(MetaTable3$PercentPOMC)*100
MetaTable3$Timepoint = x
CompilePercents = rbind(CompilePercents, MetaTable3)
}

CompilePercents = CompilePercents
CompilePercents$SampleTypeBroad = gsub(" -.*", "", CompilePercents$SampleType)
CompilePercents2 = CompilePercents %>% group_by(Timepoint, SampleTypeBroad) %>% summarise(PercentAll = sum(PercentClust_All), PercentPOMC = sum(PercentClust_POMC))

CompilePercents$SampleType = factor(CompilePercents$SampleType, levels = c("HumanFetal - CS22", "HumanFetal - GW16", "HumanFetal - GW18", "HumanFetal - GW19", "HumanFetal - GW20", "HumanFetal - GW22", "HumanFetal - GW25", "HumanAdult - Mammill", "HumanAdult - PO", "HumanAdult - SO", "HumanAdult - Tub", "MouseFetal - E10", "MouseFetal - E12", "MouseFetal - E13", "MouseFetal - E14", "MouseFetal - E15", "MouseFetal - E16", "MouseFetal - E18", "MouseAdult",    "MouseAdult - ARC",   "MouseAdult - PO",   "MouseAdult - SCN",  "MouseAdult - VMH", "MouseAdult - LH",
"MouseAdult - P4","MouseAdult - P8", "MouseAdult - P14",  "MouseAdult - P45"))
CompilePercents$Timepoint = factor(CompilePercents$Timepoint, levels = c("POMC_01", "POMC_02", "POMC_03", "POMC_04", "POMC_05", "POMC_06", "POMC_07", "POMC_08", "POMC_09", "POMC_10", "POMC_11"))
  

Fig3List[["POMC_Timing"]] = ggplot(CompilePercents, aes(fill=SampleType, y=PercentClust_All, x=Timepoint)) + geom_bar(position="stack", stat="identity") + ylab("% cells") + xlab("") + scale_y_continuous( expand= c(0,0)) + theme_classic() + theme(legend.position = "none", axis.text = element_text(size = 12), axis.title.y = element_text(size = 12), legend.title = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5)) 

############## FIGURE 3I ##############
All.integrated_POMC_Clean@meta.data$SampleTypeBroad = gsub(" -.*", "", All.integrated_POMC_Clean@meta.data$sampletype)
Idents(All.integrated_POMC_Clean) = "SampleTypeBroad"
Fig3List[["POMC_Clustering_Sample"]] = DimPlot(All.integrated_POMC_Clean)+NoLegend()

############## FIGURE 3J ##############
All.integrated_POMC_Clean@meta.data$POMC_Barcs = factor(All.integrated_POMC_Clean@meta.data$POMC_Barcs, sort(unique(All.integrated_POMC_Clean@meta.data$POMC_Barcs)))
Idents(All.integrated_POMC_Clean) = "POMC_BarcsReorder"
Fig3List[["POMC_DP"]] = DotPlot(All.integrated_POMC_Clean, rev(c("POMC","ISL1", "NDN", "SIX3","TBX3", "PCSK1", "NPY1R", "SIX6", "NKX2-1", "CRABP1", "OTP", "LEPR")), assay = "RNA")+coord_flip()+theme(axis.text.x = element_text(angle = 90,  vjust = 0.5, hjust = 0.5), legend.title = element_blank(), axis.text.y = element_text(face="italic"))+xlab("")+ylab("")

############## FIGURE 3K ##############
HumanTFs = read.csv("~/Dropbox/Columbia/Useful Gene Lists/Human TFs.csv", head=F)

TFs = intersect(rownames(All.integrated),unique(c("POMC", "AGRP", "NPY", "KISS1", "GHRH", "OXT", "AVP","CCK", HumanTFs$V1)))

HSA_NP_TF_LikelyhoodM2 = read.csv("~/Downloads/HSA_PureNeuron_NP_TF_Likelyhood_Sig_Min2.csv")
MMA_NP_TF_LikelyhoodM2 = read.csv("~/Downloads/MMA_PureNeuron_NP_TF_Likelyhood_Sig_Min2.csv")
HSA_NP_TF_LikelyhoodM2_TF = subset(HSA_NP_TF_LikelyhoodM2, HSA_NP_TF_LikelyhoodM2$X %in% TFs)
MMA_NP_TF_LikelyhoodM2_TF = subset(MMA_NP_TF_LikelyhoodM2, MMA_NP_TF_LikelyhoodM2$X %in% TFs)
HSD_NP_TF_LikelyhoodM2 = read.csv("~/Downloads/HSD_PureNeuron_NP_TF_Likelyhood_Sig_Min2.csv")
MMD_NP_TF_LikelyhoodM2 = read.csv("~/Downloads/MMD_PureNeuron_NP_TF_Likelyhood_Sig_Min2.csv")
HSD_NP_TF_LikelyhoodM2_TF = subset(HSD_NP_TF_LikelyhoodM2, HSD_NP_TF_LikelyhoodM2$X %in% TFs)
MMD_NP_TF_LikelyhoodM2_TF = subset(MMD_NP_TF_LikelyhoodM2, MMD_NP_TF_LikelyhoodM2$X %in% TFs)

PullData  <- as.data.frame(t(All.integrated@assays$RNA@counts[TFs ,]))
All.integrated.meta = All.integrated@meta.data

PullData_Norm = as.data.frame(matrix(ncol = length(PullData)+1, nrow=0))
colnames(PullData_Norm) = c(colnames(PullData), "sampletype")
for(y in unique(All.integrated.meta$sampletypeBroad)){
Pull_Subs = subset(All.integrated.meta, All.integrated.meta$sampletypeBroad == y)  
PullData_Subs = subset(PullData, row.names(PullData) %in% row.names(Pull_Subs))
PullData_Subs = apply(PullData_Subs, 2, function(x) (x-mean(x))/sd(x))
PullData_Subs = as.data.frame(PullData_Subs)
PullData_Subs$sampletype = y
PullData_Norm = rbind(PullData_Norm, PullData_Subs)
}
normalized <- function(x, ...) {(x - min(x, ...)) / (max(x, ...) - min(x, ...))}
sz =8

for(e in c(1)){
Expon = e
DoPl = list()
DoPl2 = list()
TopTFs = list()
for(g in c("POMC", "AGRP", "NPY", "KISS1", "GHRH", "OXT", "AVP","CCK")){
PullData_Norm$NewAssigns = ifelse(PullData_Norm[[g]] > 1, paste(PullData_Norm$sampletype, "+", sep=""), paste(PullData_Norm$sampletype, "-", sep=""))
NewAssignMeta = PullData_Norm %>% dplyr::select(NewAssigns)
row.names(NewAssignMeta) = row.names(PullData_Norm)
NewAssignMeta$NewAssigns = factor(NewAssignMeta$NewAssigns, rev(sort(unique(NewAssignMeta$NewAssigns))))
All.integrated = AddMetaData(All.integrated, NewAssignMeta, "NewAssignMeta")
Idents(All.integrated) = "NewAssignMeta"

HSA_NP_TF_Likelyhood_Gene = HSA_NP_TF_LikelyhoodM2_TF %>% dplyr::select(c("X", g))
colnames(HSA_NP_TF_Likelyhood_Gene) = c("X", "Gene")
HSA_NP_TF_Likelyhood_Gene = HSA_NP_TF_Likelyhood_Gene[order(-HSA_NP_TF_Likelyhood_Gene$Gene),]
HSA_NP_TF_Likelyhood_Gene$Gene[is.na(HSA_NP_TF_Likelyhood_Gene$Gene)] <- 0
HSA_NP_TF_Likelyhood_Gene$Gene = normalized(HSA_NP_TF_Likelyhood_Gene$Gene)
HSA_NP_TF_Likelyhood_Gene = HSA_NP_TF_Likelyhood_Gene[order(-HSA_NP_TF_Likelyhood_Gene$Gene),]
HSA_NP_TF_Likelyhood_Gene2 = head(HSA_NP_TF_Likelyhood_Gene, 5)

MMA_NP_TF_Likelyhood_Gene = MMA_NP_TF_LikelyhoodM2_TF %>% dplyr::select(c("X", g))
colnames(MMA_NP_TF_Likelyhood_Gene) = c("X", "Gene")
MMA_NP_TF_Likelyhood_Gene$Gene[is.na(MMA_NP_TF_Likelyhood_Gene$Gene)] <- 0
MMA_NP_TF_Likelyhood_Gene$Gene = normalized(MMA_NP_TF_Likelyhood_Gene$Gene)
MMA_NP_TF_Likelyhood_Gene = MMA_NP_TF_Likelyhood_Gene[order(-MMA_NP_TF_Likelyhood_Gene$Gene),]
MMA_NP_TF_Likelyhood_Gene2 = head(MMA_NP_TF_Likelyhood_Gene, 5)

HSD_NP_TF_Likelyhood_Gene = HSD_NP_TF_LikelyhoodM2_TF %>% dplyr::select(c("X", g))
colnames(HSD_NP_TF_Likelyhood_Gene) = c("X", "Gene")
HSD_NP_TF_Likelyhood_Gene = HSD_NP_TF_Likelyhood_Gene[order(-HSD_NP_TF_Likelyhood_Gene$Gene),]
HSD_NP_TF_Likelyhood_Gene$Gene[is.na(HSD_NP_TF_Likelyhood_Gene$Gene)] <- 0
HSD_NP_TF_Likelyhood_Gene$Gene = normalized(HSD_NP_TF_Likelyhood_Gene$Gene)
HSD_NP_TF_Likelyhood_Gene = HSD_NP_TF_Likelyhood_Gene[order(-HSD_NP_TF_Likelyhood_Gene$Gene),]
HSD_NP_TF_Likelyhood_Gene2 = head(HSD_NP_TF_Likelyhood_Gene, 5)

MMD_NP_TF_Likelyhood_Gene = MMD_NP_TF_LikelyhoodM2_TF %>% dplyr::select(c("X", g))
colnames(MMD_NP_TF_Likelyhood_Gene) = c("X", "Gene")
MMD_NP_TF_Likelyhood_Gene = MMD_NP_TF_Likelyhood_Gene[order(-MMD_NP_TF_Likelyhood_Gene$Gene),]
MMD_NP_TF_Likelyhood_Gene$Gene[is.na(MMD_NP_TF_Likelyhood_Gene$Gene)] <- 0
MMD_NP_TF_Likelyhood_Gene$Gene = normalized(MMD_NP_TF_Likelyhood_Gene$Gene)
MMD_NP_TF_Likelyhood_Gene = MMD_NP_TF_Likelyhood_Gene[order(-MMD_NP_TF_Likelyhood_Gene$Gene),]
MMD_NP_TF_Likelyhood_Gene2 = head(MMD_NP_TF_Likelyhood_Gene, 5)

GenesComb_Alt = merge(HSA_NP_TF_Likelyhood_Gene, HSD_NP_TF_Likelyhood_Gene, by = "X")
GenesComb_Alt = merge(GenesComb_Alt, MMA_NP_TF_Likelyhood_Gene, by = "X")
GenesComb_Alt = merge(GenesComb_Alt, MMD_NP_TF_Likelyhood_Gene, by = "X")
colnames(GenesComb_Alt) =c("X", "HS_Adult","HS_Fetal",  "MM_Adult", "MM_Fetal")
GenesComb_Alt$Combined = GenesComb_Alt$HS_Adult + GenesComb_Alt$HS_Fetal + GenesComb_Alt$MM_Adult + GenesComb_Alt$MM_Fetal
GenesComb_Alt = GenesComb_Alt[order(-GenesComb_Alt$Combined),]
GenesComb_Alt2 = head(GenesComb_Alt, 15)


ConservedGenes = subset(GenesComb_Alt, GenesComb_Alt$HS_Adult > mean(GenesComb_Alt$HS_Adult)+(sd(GenesComb_Alt$HS_Adult)*Expon) & GenesComb_Alt$HS_Fetal > mean(GenesComb_Alt$HS_Fetal)+(sd(GenesComb_Alt$HS_Fetal)*Expon) & GenesComb_Alt$MM_Adult > mean(GenesComb_Alt$MM_Adult)+(sd(GenesComb_Alt$MM_Adult)*Expon) & GenesComb_Alt$MM_Fetal > mean(GenesComb_Alt$MM_Fetal)+(sd(GenesComb_Alt$MM_Fetal)*Expon))

HumanOnly    = subset(GenesComb_Alt, GenesComb_Alt$HS_Adult > mean(GenesComb_Alt$HS_Adult)+(sd(GenesComb_Alt$HS_Adult)*Expon) & GenesComb_Alt$HS_Fetal > mean(GenesComb_Alt$HS_Fetal)+(sd(GenesComb_Alt$HS_Fetal)*Expon) & ! GenesComb_Alt$X %in% ConservedGenes$X)

MouseOnly    = subset(GenesComb_Alt, GenesComb_Alt$MM_Adult > mean(GenesComb_Alt$MM_Adult)+(sd(GenesComb_Alt$MM_Adult)*Expon) & GenesComb_Alt$MM_Fetal > mean(GenesComb_Alt$MM_Fetal)+(sd(GenesComb_Alt$MM_Fetal)*Expon) & ! GenesComb_Alt$X %in% ConservedGenes$X)

AdultOnly    = subset(GenesComb_Alt, GenesComb_Alt$MM_Adult > mean(GenesComb_Alt$MM_Adult)+(sd(GenesComb_Alt$MM_Adult)*Expon) & GenesComb_Alt$HS_Adult > mean(GenesComb_Alt$HS_Adult)+(sd(GenesComb_Alt$HS_Adult)*Expon) & ! GenesComb_Alt$X %in% c(ConservedGenes$X, HumanOnly$X, MouseOnly$X))

FetalOnly    = subset(GenesComb_Alt, GenesComb_Alt$MM_Fetal > mean(GenesComb_Alt$MM_Fetal)+(sd(GenesComb_Alt$MM_Fetal)*Expon) & GenesComb_Alt$HS_Fetal > mean(GenesComb_Alt$HS_Fetal)+(sd(GenesComb_Alt$HS_Fetal)*Expon) & ! GenesComb_Alt$X %in% c(ConservedGenes$X, HumanOnly$X, MouseOnly$X, AdultOnly$X))


GenesComb = c(ConservedGenes$X, HumanOnly$X, MouseOnly$X, AdultOnly$X, FetalOnly$X)

TopTFs[[paste(g, "TopCombined")]] = GenesComb

TopTFs[[paste(g, "TopEach")]] = unique(c(HSA_NP_TF_Likelyhood_Gene2$X, MMA_NP_TF_Likelyhood_Gene2$X, HSD_NP_TF_Likelyhood_Gene2$X, MMD_NP_TF_Likelyhood_Gene2$X))

MMD_Genes = subset(MMD_NP_TF_Likelyhood_Gene, MMD_NP_TF_Likelyhood_Gene$X %in% GenesComb)
MMD_Genes$sampletype = "Mouse_Fetal"
HSD_Genes = subset(HSD_NP_TF_Likelyhood_Gene, HSD_NP_TF_Likelyhood_Gene$X %in% GenesComb)
HSD_Genes$sampletype = "Human_Fetal"
MMA_Genes = subset(MMA_NP_TF_Likelyhood_Gene, MMA_NP_TF_Likelyhood_Gene$X %in% GenesComb)
MMA_Genes$sampletype = "Mouse_Adult"
HSA_Genes = subset(HSA_NP_TF_Likelyhood_Gene, HSA_NP_TF_Likelyhood_Gene$X %in% GenesComb)
HSA_Genes$sampletype = "Human_Adult"

GenesComb2 = rbind(MMD_Genes, HSD_Genes, MMA_Genes, HSA_Genes)
ForOrder = GenesComb2 %>% group_by(X) %>% dplyr::summarise(SumExp = sum(Gene))
#ForOrder = ForOrder[order(-ForOrder$SumExp),]

GenesComb2$X = factor(GenesComb2$X, levels = GenesComb)
colnames(GenesComb2) = c("Gene", "ColocScore", "sampletype")
GenesComb2$sampletype = gsub("_", "", GenesComb2$sampletype)
GenesComb2$ForMerge = paste(GenesComb2$Gene, GenesComb2$sampletype, sep="_")

PullData_Norm_Genes = PullData_Norm %>% dplyr::select(c(g, GenesComb, NewAssigns)) #Select Genes
PullData_EXPRESSION = PullData_Norm_Genes %>% group_by(NewAssigns) %>% dplyr::summarise_all(mean)
PullData_EXPRESSION = as.data.frame(PullData_EXPRESSION)
row.names(PullData_EXPRESSION) = PullData_EXPRESSION$NewAssigns
PullData_EXPRESSION$NewAssigns = NULL
PullData_EXPRESSION_gath = gather(PullData_EXPRESSION)
colnames(PullData_EXPRESSION_gath) = c("Gene", "Expression")
PullData_EXPRESSION_gath$NewAssigns = row.names(PullData_EXPRESSION) #Generate long
PullData_EXPRESSION_gath$ForMerge = paste(PullData_EXPRESSION_gath$Gene, PullData_EXPRESSION_gath$NewAssigns, sep="_")

CompileForDP2 = PullData_EXPRESSION_gath
CompileForDP2 = subset(PullData_EXPRESSION_gath, PullData_EXPRESSION_gath$NewAssigns %in% c("HumanAdult+", "HumanFetal+", "MouseAdult+", "MouseFetal+"))
CompileForDP2$NewAssigns = gsub("\\+", "", CompileForDP2$NewAssigns)
CompileForDP2$ForMerge = paste(CompileForDP2$Gene, CompileForDP2$NewAssigns, sep="_")
CompileForDP3 = merge(GenesComb2, CompileForDP2, by = "ForMerge")

CompileForDP3$Gene.x = factor(CompileForDP3$Gene.x, levels = GenesComb)

DoPl2[[g]] = ggplot(CompileForDP3, aes(x=Gene.x, y=ColocScore, color=sampletype, size = Expression)) + theme_classic() +  geom_line(aes(group = Gene.x), size=1.5, color="#F7F7F7") + geom_point(alpha = 0.7, position=position_jitter(w = 0.05, h = 0.02)) + theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size=sz, face = "italic"), plot.title = element_text(face = "bold.italic"), axis.text.y = element_text(size=sz), legend.text = element_text(size=sz), legend.title = element_blank()) +xlab("")+ylab("") + ggtitle(g)+scale_size(range = c(0, 3)) + scale_y_continuous(limits = c(0, 1.1), expand = c(0,0))

}


pdf(paste("MouseHumanPURE_ScatterONLY_Min2_SCALED_withFetal_NewOrder", e, ".pdf", sep=""), width = 10/e, height =3) 
print(DoPl2)
dev.off() 
}



## All.integrated object saved here: 'https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/Fig3_AllIntegrated_FINAL.RData'
