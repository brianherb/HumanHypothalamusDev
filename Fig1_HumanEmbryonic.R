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



EmbryonicHypo = readRDS(url('https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/HypoCS13_GW25_mt10_integrated.rds')) 

EmbryonicHypo_UMAP = as.data.frame(EmbryonicHypo@reductions$umap@cell.embeddings)

DimPlot(EmbryonicHypo)
EmbryonicHypo_Rest1 = subset(EmbryonicHypo, cells = row.names(EmbryonicHypo3D_Neurons), invert=T)
DimPlot(EmbryonicHypo_Rest1)

set.dim = 30
set.res = 1
set.kparam = 200
#ClusterFunc_All_RNA(EmbryonicHypo_Rest1)

Oligo = subset(EmbryonicHypo_UMAP, row.names(EmbryonicHypo_UMAP) %in% colnames(EmbryonicHypo_Rest1) & EmbryonicHypo_UMAP$UMAP_1 > 6 & EmbryonicHypo_UMAP$UMAP_2 < 3)
Micro = subset(EmbryonicHypo_UMAP, row.names(EmbryonicHypo_UMAP) %in% colnames(EmbryonicHypo_Rest1) & EmbryonicHypo_UMAP$UMAP_1 > 6 & EmbryonicHypo_UMAP$UMAP_2 > 3)
Endothelial = subset(EmbryonicHypo_UMAP, row.names(EmbryonicHypo_UMAP) %in% colnames(EmbryonicHypo_Rest1) & EmbryonicHypo_UMAP$UMAP_1 < -12 & EmbryonicHypo_UMAP$UMAP_2 < 3)
Blood = subset(EmbryonicHypo_UMAP, row.names(EmbryonicHypo_UMAP) %in% colnames(EmbryonicHypo_Rest1) & EmbryonicHypo_UMAP$UMAP_1 < -8 & EmbryonicHypo_UMAP$UMAP_2 > 3)
Neurons = subset(EmbryonicHypo_UMAP, row.names(EmbryonicHypo_UMAP) %in% colnames(EmbryonicHypo_Rest1) & EmbryonicHypo_UMAP$UMAP_1 > -8  & EmbryonicHypo_UMAP$UMAP_1 < 6 & EmbryonicHypo_UMAP$UMAP_2 > 3.5)

EmbryonicHypo_Rest2 = subset(EmbryonicHypo, cells = c(row.names(EmbryonicHypo3D_Neurons), row.names(Oligo), row.names(Micro), row.names(Endothelial), row.names(Blood), row.names(Neurons)), invert=T)


Endothelial_Seu = subset(EmbryonicHypo, cells = row.names(Endothelial))
DefaultAssay(Endothelial_Seu) = "RNA"
FeaturePlot(Endothelial_Seu, c("SLC38A5", "CLDN5", "ITM2A", "FLT1", "KCNJ8", "ABCC9", "ACTA2", "TAGLN"))



set.dim = 30
set.res = 0.5
set.kparam = c(50)
ClusterFunc_All_RNA(Endothelial_Seu)

DefaultAssay(Endothelial_Seu) = "integrated"
Endothelial_Seu <- FindNeighbors(Endothelial_Seu, k.param=50, dims=1:30)
Endothelial_Seu <- FindClusters(Endothelial_Seu, resolution =0.5)
Art_Endo1 = subset(Endothelial_Seu, idents= c(0)) #SLC38A5 CLDN5 ITM2A + KCNJ8, ABCC9
Art_Endo2 = subset(Endothelial_Seu, idents= c(2)) #SLC38A5 CLDN5 ITM2A + GJA4
Ven_Endo = subset(Endothelial_Seu, idents= c(1)) #SLC38A5 CLDN5 ITM2A


set.dim = 30
set.res = 1
set.kparam = 200
#ClusterFunc_All_RNA(EmbryonicHypo_Rest2)

EmbryonicHypo_Rest2 <- FindNeighbors(EmbryonicHypo_Rest2, k.param=200, dims=1:30)
EmbryonicHypo_Rest2 <- FindClusters(EmbryonicHypo_Rest2, resolution =1)
VLMC = subset(EmbryonicHypo_Rest2, idents= c(11))
VASC = subset(EmbryonicHypo_Rest2, idents= c(13))
EPENDY = subset(EmbryonicHypo_Rest2, idents= c(10))
RG = subset(EmbryonicHypo_Rest2, idents= c(0))
EmbryonicHypo_Rest3 = subset(EmbryonicHypo_Rest2, idents= c(0,10,11,13), invert=T)


set.dim = 30
set.res = 1
set.kparam = c(10,20, 50)
#ClusterFunc_All_RNA(VASC)


VASC <- FindNeighbors(VASC, k.param=50, dims=1:30)
VASC <- FindClusters(VASC, resolution =1)
vSMC = subset(VASC, idents= c(1))
Pericytes = subset(VASC, idents= c(1), invert=T)


DefaultAssay(EmbryonicHypo) = "RNA"
pdf("CheckFP_VLMC.pdf", width = 10, height = 50)
print(FeaturePlot(EmbryonicHypo, c(NewEndoList = unique(c("SLCO1C1", "SLC38A5", "MYH11", "MRC1", "CLDN5", "ITM2A", "FLT1", "CCNB1",  "CENPE", "CCL19", "MAFB", "LRRC55", "ATP1B1", "CHN2", "DDC", "LYPD1", "CDKN2B", "PLAUR", "MCAM", "PCDH17", "FBLN2", "SEMA3G", "BMX", "ABCB1", "TEK",
                 "ACTA2", "TAGLN", "ADGRF5", "EMCN", "VTN","ATP13A5", "PTH1R", "KCNJ8", "ABCC9",  "CD82", "CHST1", "PLN", "PDGFRB", "RGS5", "IGFBP2", "LUM", "DCN", "PDGFRA", "COL1A1", "COL1A2", "LUM", "DCN", "FBLN1", "TBX18", "VTN", "OGN", "IGFBP2", "IL33", "PTGDS", "NNAT", "RSPO3", "NOV", "SLC47A1", "APOD"))
)))
dev.off()


pdf("CheckFP_VLMC_v2.pdf", width = 10, height = 40)
print(FeaturePlot(EmbryonicHypo, toupper(c("Cldn5", "Adgrf5", "Emcn", "Atp13a5", "Pth1r", "Kcnj8", "Abcc9", "Apln", "Cd82", "Chst1", "Tagln", "Pln", "Bmx", "Gkn3", "Igfbp2", "Il33", "Ptgds", "Nnat", "Rspo3", "Nov", "Slc47a1", "Sox10", "Foxd3", "Aldh1a3", "Anxa11", "Vtn", "Cspg4", "Dcn", "Lum", "Pdgfra", "Slc18a2", "Klhl30", "Gfra3","Cldn19", "Mpz", "Dhrs2", "Caecam10", "Cspg5", "Olig1", "Pcdh15", "Aldoc", "Npy", "Apod"))))
dev.off()

set.dim = 30
set.res = 1
set.kparam = c(50,100)
#ClusterFunc_All_RNA(EmbryonicHypo_Rest3)

EmbryonicHypo_Rest3 <- FindNeighbors(EmbryonicHypo_Rest3, k.param=50, dims=1:30)
EmbryonicHypo_Rest3 <- FindClusters(EmbryonicHypo_Rest3, resolution =1)
EPENDY_2 = subset(EmbryonicHypo_Rest3, idents= c(16))
TANY = subset(EmbryonicHypo_Rest3, idents= c(8))
OLIGDIV = subset(EmbryonicHypo_Rest3, idents= c(13))  
BROADDIV = subset(EmbryonicHypo_Rest3, idents= c(5,6,15))
ASTROPROGEN = subset(EmbryonicHypo_Rest3, idents= c(2,10,14))
ASTRO = subset(EmbryonicHypo_Rest3, idents= c(7,12))
EGFR = subset(EmbryonicHypo_Rest3, idents= c(0))
LUM_RGS = subset(EmbryonicHypo_Rest3, idents= c(4))
EmbryonicHypo_Rest4 = subset(EmbryonicHypo_Rest3, idents= c(16,8,13,5,6,15,2,10,14,7,12,0,4), invert=T)

set.dim = 30
set.res = 1
set.kparam = c(10,50)
ClusterFunc_All_RNA(ASTRO)

ASTRO <- FindNeighbors(ASTRO, k.param=10, dims=1:30)
ASTRO <- FindClusters(ASTRO, resolution =1)
MoveToEpendy = subset(ASTRO, idents= c(2,6,7))
KeepAstro = subset(ASTRO, idents= c(2,6,7), invert=T)


set.dim = 30
set.res = 1
set.kparam = c(100)
#ClusterFunc_All_RNA(EmbryonicHypo_Rest4)

EmbryonicHypo_Rest4 <- FindNeighbors(EmbryonicHypo_Rest4, k.param=50, dims=1:30)
EmbryonicHypo_Rest4 <- FindClusters(EmbryonicHypo_Rest4, resolution =1)
LUM_RGS_2 = subset(EmbryonicHypo_Rest4, idents= c(0,3)) #LUM, RSG
NEUROPRECURSOR = subset(EmbryonicHypo_Rest4, idents= c(1,7)) 
DLX = subset(EmbryonicHypo_Rest4, idents= c(6))
NEUROPRECURSOR2 = subset(EmbryonicHypo_Rest4, idents= c(4,8)) 
OLIGPRECURSOR= subset(EmbryonicHypo_Rest4, idents= c(2,5))

#CLEAN: RadialGlia, Tanycytes, Dividing, Astrocytes, Oligodendrocyte Progenitors_1, Neural Precursor_1, Oligodendrocyte Progenitors_2
MoveRGtoOP2 = subset(EmbryonicHypo_UMAP, row.names(EmbryonicHypo_UMAP) %in% colnames(RG) & EmbryonicHypo_UMAP$UMAP_1 < 0 & EmbryonicHypo_UMAP$UMAP_1 > -5 )
MoveRGtoDiv = subset(EmbryonicHypo_UMAP, row.names(EmbryonicHypo_UMAP) %in% colnames(RG)  & EmbryonicHypo_UMAP$UMAP_1 < -5 )
RG_Keep = subset(RG, cells = c(row.names(MoveRGtoOP2), row.names(MoveRGtoDiv)), invert=T)

MoveTanytoOP1 = subset(EmbryonicHypo_UMAP, row.names(EmbryonicHypo_UMAP) %in% colnames(TANY) & EmbryonicHypo_UMAP$UMAP_2 > -3)
MoveTanytoOP2 = subset(EmbryonicHypo_UMAP, row.names(EmbryonicHypo_UMAP) %in% colnames(TANY) & EmbryonicHypo_UMAP$UMAP_1 < -2 & EmbryonicHypo_UMAP$UMAP_2 < 3)
Tany_Keep = subset(TANY, cells = c(row.names(MoveTanytoOP1), row.names(MoveTanytoOP2)), invert=T)

ReclusterDiv = subset(EmbryonicHypo_UMAP, row.names(EmbryonicHypo_UMAP) %in% colnames(BROADDIV) & EmbryonicHypo_UMAP$UMAP_1 > -3)
BroadDiv_Keep = subset(BROADDIV, cells = c(row.names(ReclusterDiv)), invert=T)

MoveAstrotoEpendy = subset(EmbryonicHypo_UMAP, row.names(EmbryonicHypo_UMAP) %in% colnames(KeepAstro) & EmbryonicHypo_UMAP$UMAP_2 < -9.5)
MoveAstrotoAstProgen = subset(EmbryonicHypo_UMAP, row.names(EmbryonicHypo_UMAP) %in% colnames(KeepAstro) & EmbryonicHypo_UMAP$UMAP_2 > -6 & EmbryonicHypo_UMAP$UMAP_2 < -4)
MoveAstrotoRG = subset(EmbryonicHypo_UMAP, row.names(EmbryonicHypo_UMAP) %in% colnames(KeepAstro) & EmbryonicHypo_UMAP$UMAP_2 > -4)
Astro_KeepV2 = subset(KeepAstro, cells = c(row.names(MoveAstrotoEpendy), row.names(MoveAstrotoAstProgen), row.names(MoveAstrotoRG)), invert=T)

MoveOP1toRG = subset(EmbryonicHypo_UMAP, row.names(EmbryonicHypo_UMAP) %in% colnames(EGFR) & EmbryonicHypo_UMAP$UMAP_1 > 3 & EmbryonicHypo_UMAP$UMAP_2 > -5)
MoveOP1toEpendy = subset(EmbryonicHypo_UMAP, row.names(EmbryonicHypo_UMAP) %in% colnames(EGFR) & EmbryonicHypo_UMAP$UMAP_1 > 3 & EmbryonicHypo_UMAP$UMAP_2 < -5)
MoveOP1toAstProgen = subset(EmbryonicHypo_UMAP, row.names(EmbryonicHypo_UMAP) %in% colnames(EGFR) & EmbryonicHypo_UMAP$UMAP_1 > 2 & EmbryonicHypo_UMAP$UMAP_2 < -5 & EmbryonicHypo_UMAP$UMAP_2 > -10)
MoveOP1toNeurProgen = subset(EmbryonicHypo_UMAP, row.names(EmbryonicHypo_UMAP) %in% colnames(EGFR) & EmbryonicHypo_UMAP$UMAP_1 < 2 & EmbryonicHypo_UMAP$UMAP_2 < -4)
OP1_Keep = subset(EGFR, cells = c(row.names(MoveOP1toRG), row.names(MoveOP1toAstProgen), row.names(MoveOP1toNeurProgen), row.names(MoveOP1toEpendy)), invert=T)

MoveNP1toAstProgen = subset(EmbryonicHypo_UMAP, row.names(EmbryonicHypo_UMAP) %in% colnames(NEUROPRECURSOR) & EmbryonicHypo_UMAP$UMAP_2 < -7)
MoveNP1toRG = subset(EmbryonicHypo_UMAP, row.names(EmbryonicHypo_UMAP) %in% colnames(NEUROPRECURSOR) & EmbryonicHypo_UMAP$UMAP_1 > 2 & EmbryonicHypo_UMAP$UMAP_2 > -5)
MoveNP1toOP1 = subset(EmbryonicHypo_UMAP, row.names(EmbryonicHypo_UMAP) %in% colnames(NEUROPRECURSOR) & EmbryonicHypo_UMAP$UMAP_2 > 0)
MoveNP1toOP2 = subset(EmbryonicHypo_UMAP, row.names(EmbryonicHypo_UMAP) %in% colnames(NEUROPRECURSOR) & EmbryonicHypo_UMAP$UMAP_1 < -2)
NP1_Keep = subset(NEUROPRECURSOR, cells = c(row.names(MoveNP1toAstProgen), row.names(MoveNP1toRG), row.names(MoveNP1toOP1), row.names(MoveNP1toOP2)), invert=T)

MoveOP2toAstProgen = subset(EmbryonicHypo_UMAP, row.names(EmbryonicHypo_UMAP) %in% colnames(OLIGPRECURSOR) & EmbryonicHypo_UMAP$UMAP_2 < -7)
MoveOP2toRG = subset(EmbryonicHypo_UMAP, row.names(EmbryonicHypo_UMAP) %in% colnames(OLIGPRECURSOR) & EmbryonicHypo_UMAP$UMAP_1 > 0.5 & EmbryonicHypo_UMAP$UMAP_2 > -5)
OP2_Keep = subset(OLIGPRECURSOR, cells = c(row.names(MoveOP2toAstProgen), row.names(MoveOP2toRG)), invert=T)

ReclusterDivSeu = subset(BROADDIV, cells = row.names(ReclusterDiv))
set.dim = 30
set.res = 1
set.kparam = c(20, 50)
ClusterFunc_All_RNA(ReclusterDivSeu)

ReclusterDivSeu <- FindNeighbors(ReclusterDivSeu, k.param=20, dims=1:30)
ReclusterDivSeu <- FindClusters(ReclusterDivSeu, resolution =1)
MovetoOP2 = subset(ReclusterDivSeu, idents= c(0,3))
MovetoNP1 = subset(ReclusterDivSeu, idents= c(1,2))
MovetoOP1 = subset(ReclusterDivSeu, idents= c(4))

#Clean Neurons
CombEmbryonicNeur = subset(EmbryonicHypo, cells = c(row.names(EmbryonicHypo3D_Neurons), row.names(Neurons)))
CombEmbryonicNeur_clean = subset(EmbryonicHypo_UMAP, row.names(EmbryonicHypo_UMAP) %in% colnames(CombEmbryonicNeur) & EmbryonicHypo_UMAP$UMAP_1 > -7.5 & EmbryonicHypo_UMAP$UMAP_1 < 5 & EmbryonicHypo_UMAP$UMAP_2 > -2.5)
CleanedNeurons = subset(EmbryonicHypo, cells = c(row.names(CombEmbryonicNeur_clean)))

ToReassign = subset(CombEmbryonicNeur, cells = colnames(CleanedNeurons), invert=T)
ExtraBlood = subset(EmbryonicHypo_UMAP, row.names(EmbryonicHypo_UMAP) %in% colnames(ToReassign) & EmbryonicHypo_UMAP$UMAP_1 < -5 & EmbryonicHypo_UMAP$UMAP_2 > 5)
ExtraArt1 = subset(EmbryonicHypo_UMAP, row.names(EmbryonicHypo_UMAP) %in% colnames(ToReassign) & EmbryonicHypo_UMAP$UMAP_1 < -5 & EmbryonicHypo_UMAP$UMAP_2 < 5 & EmbryonicHypo_UMAP$UMAP_2 > -4)
ExtraOlig = subset(EmbryonicHypo_UMAP, row.names(EmbryonicHypo_UMAP) %in% colnames(ToReassign) & EmbryonicHypo_UMAP$UMAP_1 > 5)
ExtraDivOlig = subset(EmbryonicHypo_UMAP, row.names(EmbryonicHypo_UMAP) %in% colnames(ToReassign) & EmbryonicHypo_UMAP$UMAP_1 < -5 &  EmbryonicHypo_UMAP$UMAP_2 < -4)
ExtraEpendy = subset(EmbryonicHypo_UMAP, row.names(EmbryonicHypo_UMAP) %in% colnames(ToReassign) & EmbryonicHypo_UMAP$UMAP_2 < -7)
ExtraAstroProgen = subset(EmbryonicHypo_UMAP, row.names(EmbryonicHypo_UMAP) %in% colnames(ToReassign) & EmbryonicHypo_UMAP$UMAP_2 > -7 & EmbryonicHypo_UMAP$UMAP_1 > -5 & EmbryonicHypo_UMAP$UMAP_1 < 8)



##Cleaned
RG_Clean = subset(EmbryonicHypo, cells = c(colnames(RG_Keep), row.names(MoveAstrotoRG), row.names(MoveOP1toRG), row.names(MoveNP1toRG), row.names(MoveOP2toRG)))
Tany_Clean = subset(EmbryonicHypo, cells =  c(colnames(Tany_Keep)))
BroadDiv_Clean = subset(EmbryonicHypo, cells =  c(colnames(BroadDiv_Keep), row.names(MoveRGtoDiv)))
OP1_Clean = subset(EmbryonicHypo, cells =  c(colnames(OP1_Keep), row.names(MoveTanytoOP1), row.names(MoveNP1toOP1), colnames(MovetoOP1)))
OP2_Clean = subset(EmbryonicHypo, cells =  c(colnames(OP2_Keep), row.names(MoveRGtoOP2),  row.names(MoveTanytoOP2), row.names(MoveNP1toOP2), colnames(MovetoOP2)))
NP1_Clean = subset(EmbryonicHypo, cells =  c(colnames(NP1_Keep), row.names(MoveOP1toNeurProgen), colnames(MovetoNP1)))
Astro_Clean = subset(EmbryonicHypo, cells =  c(colnames(Astro_KeepV2)))
Ependy_Clean = subset(EmbryonicHypo, cells =  c(row.names(MoveAstrotoEpendy), row.names(MoveOP1toEpendy), colnames(EPENDY), colnames(EPENDY_2), colnames(MoveToEpendy), row.names(ExtraEpendy)))
AstroProgen_Clean = subset(EmbryonicHypo, cells =  c(row.names(MoveAstrotoAstProgen), row.names(MoveOP1toAstProgen), row.names(MoveNP1toAstProgen), row.names(MoveOP2toAstProgen), colnames(ASTROPROGEN), row.names(ExtraAstroProgen)))
Blood_Clean = subset(EmbryonicHypo, cells = c(row.names(Blood), row.names(ExtraBlood)))
Art1_Clean = subset(EmbryonicHypo, cells = c(colnames(Art_Endo1), row.names(ExtraArt1)))
Oligodendrocytes_Clean = subset(EmbryonicHypo, cells = c(row.names(ExtraOlig), row.names(Oligo)))
DivOlig_Clean = subset(EmbryonicHypo, cells = c(colnames(OLIGDIV), row.names(ExtraDivOlig)))


#Cluster Oligodendrocytes
set.dim = 30
set.res = 1
set.kparam = c(20, 50, 100)
ClusterFunc_All_RNA(Oligodendrocytes_Clean)

Oligodendrocytes_Clean <- FindNeighbors(Oligodendrocytes_Clean, k.param=50, dims=1:30)
Oligodendrocytes_Clean <- FindClusters(Oligodendrocytes_Clean, resolution =1)
Olig_Mat = subset(Oligodendrocytes_Clean, idents= c(4))
Olig_Maturing = subset(Oligodendrocytes_Clean, idents= c(3,6))
Olig_Immat = subset(Oligodendrocytes_Clean, idents= c(4, 3, 6), invert=T)



### COMPILE ASSIGNMENTS
MetaDatasets = list( Micro, VLMC, Art1_Clean, Art_Endo2, Ven_Endo, Blood_Clean, CleanedNeurons, vSMC, Pericytes, Ependy_Clean, RG_Clean, Tany_Clean, DivOlig_Clean, BroadDiv_Clean, AstroProgen_Clean, Astro_Clean, OP1_Clean, LUM_RGS, LUM_RGS_2, NP1_Clean, DLX, NEUROPRECURSOR2, OP2_Clean, Olig_Mat, Olig_Maturing, Olig_Immat)
names(MetaDatasets)  = c( "Microglia", "Pericytes_1", "Endothelial [Arterial_1]", "Endothelial [Arterial_2]", "Endothelial [Venous]", "Blood", "Neurons", "vSMC", "Pericytes_2", "Ependymal", "RadialGlia", "Tanycytes", "Oligodendrocytes [Dividing]", "Dividing", "Astrocyte Progenitors", "Astrocytes", "Oligodendrocyte Progenitors_1", "VLMC", "VLMC", "Neural Progenitors_1", "Neural Progenitors_2", "Neural Progenitors_2", "Oligodendrocyte Progenitors_2", "Oligodendrocytes [Mature]", "Oligodendrocytes [Maturing]", "Oligodendrocytes [Immature]")
CurrentMeta = GenerateMetaData(MetaDatasets)



### FIGURE 1

EmbryonicHypo@meta.data$Timepoint = gsub("_.*", "", gsub("T", "", EmbryonicHypo@meta.data$sample))
EmbryonicHypo3D$Timepoint = gsub("_.*", "", gsub("T", "", row.names(EmbryonicHypo3D)))
EmbryonicHypo3D$Barcodes = gsub("-.*", "", gsub("T_hypo1", "", gsub("3V_", "",  gsub("_34", "",  gsub("_2_", "_", row.names(EmbryonicHypo3D))))))
EmbryonicHypo3D$Barcodes = gsub("_hypo", "", EmbryonicHypo3D$Barcodes)


############## FIGURE 1B ##############
Idents(EmbryonicHypo) = "Timepoint"
p = plot_ly(EmbryonicHypo3D, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, size = 1, color = ~Timepoint, colors = "Spectral") ##PLOT IN PAPER - FIG 1B
htmlwidgets::saveWidget(p, "./Figure1_3D_BySample_1APR22.html")

############## SUPP FIG 1 ##############
pdf("Fig1Supp_SamplesSplit.pdf", width = 10, height = 12)
print(DimPlot(EmbryonicHypo, split.by = "Timepoint", ncol = 3, cols = c("#9e0142", "#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#e6f598", "#abdda4", "#66c2a5", "#3288bd", "#5e4fa2")))
dev.off()


ColPal = c("RadialGlia" = "#9e1b44", "Dividing" = "#d64050", "Neurons" = "#f36c44", "Neural Progenitors_1" = "#faad60", "Neural Progenitors_2" = "#fee08b", "Astrocyte Progenitors" = "#e4eb9a", "Astrocytes" = "#acd7a5", "Ependymal" = "#7eb996", "Tanycytes" = "#537a63", "Oligodendrocytes [Mature]" = "#225b7e", "Oligodendrocytes [Maturing]" = "#2d75a4", "Oligodendrocytes [Immature]" = "#42a0d9", "OligodendrocyteProgenitors_1" = "#5ebadf", "Oligodendrocyte Progenitors_2" = "#7dcbea", "Oligodendrocytes [Dividing]" = "#81afdd", "Pericytes_1" = "#8d7cba", "Pericytes_2" = "#6e5aa7", "Endothelial [Arterial_1]" = "#7d4a9d", "Endothelial [Arterial_2]" = "#9d5ba5", "Endothelial [Venous]" = "#b781b9", "vSMC" = "#ecbed8", "VLMC" = "#cc6cab", "Microglia" = "#b8258f", "Blood" = "#731e5f")


############## FIGURE 1C ############## 
Assignw3D = merge(EmbryonicHypo3D, CurrentMeta, by = 0)
Assignw3D$Pop = factor(Assignw3D$Pop, levels = c("RadialGlia", "Dividing", "Neurons", "Neural Progenitors_1", "Neural Progenitors_2", "Astrocyte Progenitors", "Astrocytes", "Ependymal", "Tanycytes", "Oligodendrocytes [Mature]", "Oligodendrocytes [Maturing]", "Oligodendrocytes [Immature]", "OligodendrocyteProgenitors_1", "Oligodendrocyte Progenitors_2", "Oligodendrocytes [Dividing]", "Pericytes_1", "Pericytes_2", "Endothelial [Arterial_1]", "Endothelial [Arterial_2]", "Endothelial [Venous]", "vSMC", "VLMC", "Microglia", "Blood"))
Assignw3D = Assignw3D[order(Assignw3D$Pop), ]

p = plot_ly(Assignw3D, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, size = 1, color = ~Pop, colors = ColPal, type="scatter3d")
htmlwidgets::saveWidget(p, "./Figure1_3D_AllAssignments_2MAY22.html")

############## FIGURE 1D ############## 
EmbryonicHypoMeta = EmbryonicHypo@meta.data

CompilePercents = as.data.frame(matrix(ncol =3, nrow=0))
colnames(CompilePercents) = c("Var1", "Percent", "Timepoint")
for(x in unique(EmbryonicHypoMeta$Timepoint)){
Subs= subset(EmbryonicHypoMeta, EmbryonicHypoMeta$Timepoint == x) 
MetaTable = as.data.frame(table(Subs$CurrentMeta))
MetaTable$Percent = MetaTable$Freq/dim(Subs)[1]*100
MetaTable$Freq = NULL
MetaTable$Timepoint = x
CompilePercents = rbind(CompilePercents, MetaTable)
}
CompilePercents2 = merge(CompilePercents, ColPal, by.x = "Var1", by.y = 0)
CompilePercents2$Var1 = factor(CompilePercents2$Var1, levels = rev(c("RadialGlia", "Dividing", "Neurons", "Neural Progenitors_1", "Neural Progenitors_2", "Astrocyte Progenitors", "Astrocytes", "Ependymal", "Tanycytes", "Oligodendrocytes [Mature]", "Oligodendrocytes [Maturing]", "Oligodendrocytes [Immature]", "OligodendrocyteProgenitors_1", "Oligodendrocyte Progenitors_2", "Oligodendrocytes [Dividing]", "Pericytes_1", "Pericytes_2", "Endothelial [Arterial_1]", "Endothelial [Arterial_2]", "Endothelial [Venous]", "vSMC", "VLMC", "Microglia", "Blood")))
CompilePercents2$Timepoint = factor(CompilePercents2$Timepoint, levels = c("CS13", "CS14", "CS15","CS22", "GW16",  "GW18", "GW19", "GW20",  "GW22", "GW25"))
  
StackedPlot = ggplot(CompilePercents2, aes(fill=Var1, y=Percent, x=Timepoint)) + geom_bar(position="stack", stat="identity") + ylab("% cells") + xlab("") + scale_y_continuous(limits = c(0,100), expand= c(0,0)) + theme_classic() + theme(legend.position = "right", axis.text = element_text(size = 16), axis.title.y = element_text(size = 16), legend.title = element_blank()) + scale_fill_manual(values = ColPal)

pdf("Fetal_ForPaper_StackedbyPercent.pdf", width = 18*1.3, height = 3*1.3)
print(StackedPlot)
dev.off()


############## FIGURE 1E ############## 
MetaDatasets_Neuro = list( Micro, VLMC, Art1_Clean, Art_Endo2, Ven_Endo, Blood_Clean, CleanedNeurons, vSMC, Pericytes, Ependy_Clean, RG_Clean, Tany_Clean, DivOlig_Clean, BroadDiv_Clean, AstroProgen_Clean, Astro_Clean, OP1_Clean, LUM_RGS, LUM_RGS_2, NP1_Clean, DLX, NEUROPRECURSOR2, OP2_Clean, Olig_Mat, Olig_Maturing, Olig_Immat)
names(MetaDatasets_Neuro)  = c("Other", "Other", "Other", "Other", "Other", "Other", "Neurons", "Other", "Other", "Other", "RadialGlia", "Other", "Other", "Dividing", "Other", "Other", "Other", "Other", "Other", "Neural Progenitors_1", "Neural Progenitors_2", "Neural Progenitors_2", "Other", "Other", "Other", "Other")
CurrentMeta_Neuro = GenerateMetaData(MetaDatasets_Neuro)
ColPal_Neuro = c("RadialGlia" = "#9e1b44", "Dividing" = "#b8258f", "Neural Progenitors_1" = "#7d4a9d", "Neural Progenitors_2" = "#7dcbea", "Neurons" = "#7eb996",  "Other" = "lightGrey")
CurrentMeta_Neuro$Pop = factor(CurrentMeta_Neuro$Pop, levels = c("RadialGlia", "Dividing", "Neural Progenitors_1", "Neural Progenitors_2", "Neurons"))
Assignw3D_Neuro = merge(EmbryonicHypo3D, CurrentMeta_Neuro, by = 0)
p = plot_ly(Assignw3D_Neuro, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, size = 1, color = ~Pop, colors = ColPal_Neuro)
htmlwidgets::saveWidget(p, "./Figure1_3D_Neuro_Trajectory_2MAY22.html")

############## FIGURE 1F ############## 
MetaDatasets_Astro = list( Micro, VLMC, Art1_Clean, Art_Endo2, Ven_Endo, Blood_Clean, CleanedNeurons, vSMC, Pericytes, Ependy_Clean, RG_Clean, Tany_Clean, DivOlig_Clean, BroadDiv_Clean, AstroProgen_Clean, Astro_Clean, OP1_Clean, LUM_RGS, LUM_RGS_2, NP1_Clean, DLX, NEUROPRECURSOR2, OP2_Clean, Olig_Mat, Olig_Maturing, Olig_Immat)
names(MetaDatasets_Astro)  = c("Other", "Other", "Other", "Other", "Other", "Other", "Other", "Other", "Other", "Ependymal", "RadialGlia", "Tanycytes", "Other", "Dividing", "Astrocyte Progenitors", "Astrocytes", "Other", "Other", "Other", "Other", "Other", "Other", "Other", "Other", "Other", "Other")
CurrentMeta_Astro = GenerateMetaData(MetaDatasets_Astro)
ColPal_Astro = c("RadialGlia" = "#9e1b44", "Dividing" = "#b8258f", "Astrocyte Progenitors" = "#7d4a9d", "Astrocytes" = "#7dcbea", "Tanycytes" = "#2d75a4",  "Ependymal" = "#7eb996")
CurrentMeta_Astro$Pop = factor(CurrentMeta_Astro$Pop, levels = c("RadialGlia",  "Dividing", "Astrocyte Progenitors", "Astrocytes", "Tanycytes",  "Ependymal"))
Assignw3D_Astro = merge(EmbryonicHypo3D, CurrentMeta_Astro, by = 0)
p = plot_ly(Assignw3D_Astro, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, size = 1, color = ~Pop, colors = ColPal_Astro)
htmlwidgets::saveWidget(p, "./Figure1_3D_Astro_Trajectory_2MAY22.html")


############## FIGURE 1G ############## 
MetaDatasets_Olig = list( Micro, VLMC, Art1_Clean, Art_Endo2, Ven_Endo, Blood_Clean, CleanedNeurons, vSMC, Pericytes, Ependy_Clean, RG_Clean, Tany_Clean, DivOlig_Clean, BroadDiv_Clean, AstroProgen_Clean, Astro_Clean, OP1_Clean, LUM_RGS, LUM_RGS_2, NP1_Clean, DLX, NEUROPRECURSOR2, OP2_Clean, Olig_Mat, Olig_Maturing, Olig_Immat)
names(MetaDatasets_Olig)  = c("Other", "Other", "Other", "Other", "Other", "Other", "Neurons", "Other", "Other", "Other", "RadialGlia", "Other", "Oligodendrocytes [Dividing]", "Dividing", "Other", "Other", "Oligodendrocyte Progenitors_1", "Other", "Other", "Neural Progenitors_1", "Other", "Other", "Oligodendrocyte Progenitors_2", "Oligodendrocytes [Mature]", "Oligodendrocytes [Maturing]", "Oligodendrocytes [Immature]")
CurrentMeta_Olig = GenerateMetaData(MetaDatasets_Olig)
ColPal_Olig = c("RadialGlia" = "#9e1b44", "Dividing" = "#b8258f", "Oligodendrocytes [Dividing]" = "#7d4a9d", "Oligodendrocyte Progenitors_1" = "#7dcbea", "Oligodendrocyte Progenitors_2" =  "#2d75a4", "Oligodendrocytes [Immature]" = "#7eb996",  "Oligodendrocytes [Maturing]" = "#fee08b",  "Oligodendrocytes [Mature]" = "#f36c44")
                
CurrentMeta_Olig$Pop = factor(CurrentMeta_Olig$Pop, levels = c("RadialGlia",  "Dividing", "Oligodendrocytes [Dividing]","Oligodendrocyte Progenitors_1", "Oligodendrocyte Progenitors_2","Oligodendrocytes [Immature]","Oligodendrocytes [Maturing]",  "Oligodendrocytes [Mature]", "Other"))
Assignw3D_Olig = merge(EmbryonicHypo3D, CurrentMeta_Olig, by = 0)
p = plot_ly(Assignw3D_Olig, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, size = 1, color = ~Pop, colors = ColPal_Olig)
htmlwidgets::saveWidget(p, "./Figure1_3D_Olig_Trajectory_2MAY22.html")



############## FIGURE 1H ############## 
CurrentMeta$Pop = factor(CurrentMeta$Pop, levels = rev(c("RadialGlia", "Dividing", "Neurons", "Neural Progenitors_1", "Neural Progenitors_2", "Astrocyte Progenitors", "Astrocytes", "Ependymal", "Tanycytes", "Oligodendrocytes [Mature]", "Oligodendrocytes [Maturing]", "Oligodendrocytes [Immature]", "OligodendrocyteProgenitors_1", "Oligodendrocyte Progenitors_2", "Oligodendrocytes [Dividing]", "Pericytes_1", "Pericytes_2", "Endothelial [Arterial_1]", "Endothelial [Arterial_2]", "Endothelial [Venous]", "vSMC", "VLMC", "Microglia", "Blood")))
EmbryonicHypo = AddMetaData(EmbryonicHypo, CurrentMeta, "CurrentMeta")
Idents(EmbryonicHypo) = "CurrentMeta"

DP = Seurat::DotPlot(EmbryonicHypo, assay = "RNA", features= c("HOPX", "MKI67", "STMN2", "SYT1", "PAX6", "DLX5", "NKX2−1", "GFAP", "AQP4", "SLC1A3", "FOXJ1", "CCDC153", "RFX3", "RFX4", "CRYM", "FYN", "MBP", "CNP", "TNS3", "TCF7L2", "PLP1", "MAG", "MOG", "PDGFRA", "OLIG2", "APOD", "OLIG1", "EGFR", "VCAN", "SLC7A11", "COL1A1", "COL3A1", "DCN", "RGS5", "PDGFRB", "KCNJ8", "ABCC9", "CLDN5", "ITM2A", "FLT1", "ACTA2", "TAGLN", "LUM", "C1QB", "AIF1", "HBA2"),  cols = c("lightgrey", "navy")) + theme(axis.text.x = element_text(size=10, face="italic", angle = 90, vjust = 0.5, hjust=1), axis.text.y = element_text(size=10), legend.title = element_text(size=8)) + ylab("")

pdf("Fetal_ForPaper_Fig1DotPlot.pdf", width = 16, height = 5)
print(DP)
dev.off()

DP = Seurat::DotPlot(EmbryonicHypo, assay = "RNA", features= SHH_Go,  cols = c("#eeefef", "#2d75a4")) + theme(axis.text.x = element_text(size=10, face="italic", angle = 90, vjust = 0.5, hjust=1), axis.text.y = element_text(size=10), legend.title = element_text(size=8)) + ylab("")

pdf("Fetal_ForPaper_Fig1_SHHGenes.pdf", width = 16, height = 5)
print(DP)
dev.off()

## EmbryonicHypo object saved here: 'https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/Fig1_EmbryonicHypo_FINAL.RData'