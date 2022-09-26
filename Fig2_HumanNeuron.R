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


### FIGURE 2
load(url('https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/EdKaHypoNeurons_mt10_3D_integrated.rda')) # EdKaHypoNeurons_mt10_3D


DefaultAssay(EdKaHypoNeurons_mt10_3D) = "integrated"
EdKaHypoNeurons_mt10_3D <- FindNeighbors(EdKaHypoNeurons_mt10_3D, k.param=100, dims=1:30)
EdKaHypoNeurons_mt10_3D <- FindClusters(EdKaHypoNeurons_mt10_3D, resolution =1)
Nac_Cells = subset(EdKaHypoNeurons_mt10_3D, idents= c(6,13,24,28))
AdultNeuroExNac = subset(EdKaHypoNeurons_mt10_3D, idents= c(6,13,24,28), invert=T)

EdKaHypoNeurons_mt10_3D_UMAP = as.data.frame(EdKaHypoNeurons_mt10_3D@reductions$umap@cell.embeddings)
Nac_Cells_Pt2 = subset(EdKaHypoNeurons_mt10_3D_UMAP, row.names(EdKaHypoNeurons_mt10_3D_UMAP) %in% colnames(AdultNeuroExNac) & EdKaHypoNeurons_mt10_3D_UMAP$UMAP_1 > -5 & EdKaHypoNeurons_mt10_3D_UMAP$UMAP_1 < -0.2 & EdKaHypoNeurons_mt10_3D_UMAP$UMAP_2 > 1 & EdKaHypoNeurons_mt10_3D_UMAP$UMAP_2 < 6 | row.names(EdKaHypoNeurons_mt10_3D_UMAP) %in% colnames(AdultNeuroExNac) & EdKaHypoNeurons_mt10_3D_UMAP$UMAP_1 > -5 & EdKaHypoNeurons_mt10_3D_UMAP$UMAP_1 < 0.3 & EdKaHypoNeurons_mt10_3D_UMAP$UMAP_2 > 3 & EdKaHypoNeurons_mt10_3D_UMAP$UMAP_2 < 5.5) 
AdultNeuroExNac_V2 = subset(AdultNeuroExNac, cells = row.names(Nac_Cells_Pt2), invert=T)


Idents(AdultNeuroExNac_V2) = "hicat_clust"
Hicat = as.data.frame(table(AdultNeuroExNac_V2@meta.data$hicat_clust))
write.csv(Hicat, "Hicat_Clusts.csv")
Hicat_Assigns = read.csv("~/Downloads/Hicat_Assigns.csv")
Hicats_Nuclei = as.data.frame(matrix(ncol = 2, nrow=0))
colnames(Hicats_Nuclei) = c("Hicat_Nuclei", "Hicat_Subtype")
Cumulative = c()
for(x in unique(Hicat_Assigns$X)){
PullHicat = subset(AdultNeuroExNac_V2@meta.data, AdultNeuroExNac_V2@meta.data$hicat_clust %in% x)
Pull2 = subset(Hicat_Assigns, Hicat_Assigns$X %in% x)
Cumulative = c(Cumulative, Pull2$Assignments)
CumulativeCount = subset(Cumulative, Cumulative == Pull2$Assignments)
Outs = as.data.frame(rep(paste(Pull2$Assignments, length(CumulativeCount), sep="_"), dim(PullHicat)[1]))
colnames(Outs) = "Hicat_Nuclei"
Outs$Hicat_Subtype = Pull2$Subtype.Exp
row.names(Outs) = row.names(PullHicat)
Hicats_Nuclei = rbind(Hicats_Nuclei, Outs)
}

Hicat_Nuclei2 = Hicats_Nuclei
Hicat_Nuclei2$Hicat_Subtype = NULL
AdultNeuroExNac_V2 = AddMetaData(AdultNeuroExNac_V2, Hicat_Nuclei2, "Hicat_Nuclei2")


############## FIGURE 2D ############## 
Hicat_Nuclei3 = subset(Hicat_Nuclei2, ! Hicat_Nuclei2$Hicat_Nuclei %in% c("Ukn_1", "Ukn_10", "Ukn_11", "Ukn_12", "Ukn_2", "Ukn_3", "Ukn_4", "Ukn_5", "Ukn_6", "Ukn_7","Ukn_8", "Ukn_9", "Excl (< 10 cells)_1",  "Excl (< 10 cells)_10", "Excl (< 10 cells)_11", "Excl (< 10 cells)_2",  "Excl (< 10 cells)_3", "Excl (< 10 cells)_4",  "Excl (< 10 cells)_5",  "Excl (< 10 cells)_6",  "Excl (< 10 cells)_7",  "Excl (< 10 cells)_8",  "Excl (< 10 cells)_9"))
HiCatSeu =subset(AdultNeuroExNac_V2, cells = row.names(Hicat_Nuclei3))
HicatData = as.data.frame(t(HiCatSeu@assays$RNA@counts))
HicatData2 = HicatData %>% dplyr::select("SLC32A1", "LHX8", "PMCH", "PNOC", "NR2F2", "NR2F1", "VIP", "SIX3", "SIX6", "POU3F2", "AVP", "CRH", "OXT", "SIM1", "LHX9", "IRX5", "FOXA1", "BARHL1", "LMX1A", "LHX6", "ARX", "MEIS2", "PITX2", "LHX1", "FOXB1", "GHRH", "SOX4", "HMX2", "LEPR", "NPY", "AGRP", "POMC", "TBX3", "HDC", "NR5A1", "NPTX2", "FEZF1", "SLC17A6", "GAD2", "HCRT", "CARTPT", "RELN", "HMX3", "VAX1", "GABRA5", "SOX14", "GSX1", "GAD1", "FEZF1", "TBR1")
HicatData2 = apply(HicatData2, 2, scale)
HicatData2 = as.data.frame(HicatData2)
row.names(HicatData2) = row.names(HicatData)
HicatData3 = merge(HicatData2, Hicat_Nuclei3, by = 0)
HicatData3$Row.names = NULL
HicatData4 = HicatData3 %>% group_by(Hicat_Nuclei) %>% dplyr::summarise(across(everything(), mean))
HicatData4 = as.data.frame(HicatData4)
row.names(HicatData4) = HicatData4$Hicat_Nuclei
HicatData4$Hicat_Nuclei = NULL
Hicat_gath = gather(HicatData4, Gene, Expression)
Hicat_gath$Cluster = row.names(HicatData4)
Hicat_gath = subset(Hicat_gath, Hicat_gath$Gene %in% c("MEIS2", "AGRP", "POMC", "TBX3", "LHX8", "ARX", "LHX6", "LHX9", "FOXB1", "PITX2", "SIM1", "POU3F2", "LMX1A", "VIP", "HDC", "NR5A1", "FEZF1", "TBR1", "SOX14"))
Hicat_gath$Log2 = log2(Hicat_gath$Expression+1)
ExtraData = as.data.frame(matrix(ncol = 0, nrow = 12*19))
ExtraData$Gene = c("MEIS2", "AGRP", "POMC", "TBX3", "LHX8", "ARX", "LHX6", "LHX9", "FOXB1", "PITX2", "SIM1", "POU3F2", "LMX1A", "VIP", "HDC", "NR5A1", "FEZF1", "TBR1", "SOX14")
ExtraData$Expression = NA
ExtraData$Log2 = NA
ExtraData$Cluster = c(rep("NA1", 19), rep("NA2", 19), rep("NA3", 19), rep("NA4", 19), rep("NA5", 19), rep("NA6", 19), rep("NA7", 19), rep("NA8", 19), rep("NA9", 19), rep("NA10", 19), rep("NA11", 19), rep("NA12", 19))
Hicat_gath2 = rbind(Hicat_gath, ExtraData)

Hicat_gath2$Gene = factor(Hicat_gath2$Gene, levels = rev(c("MEIS2", "AGRP", "POMC", "TBX3", "LHX8", "ARX", "LHX6", "LHX9", "FOXB1", "PITX2", "SIM1", "POU3F2", "LMX1A", "VIP", "HDC", "NR5A1", "FEZF1", "TBR1", "SOX14")))
Hicat_gath2$Cluster = factor(Hicat_gath2$Cluster, levels = c("AH_1", "AH_2", "AH_3", "AH_4", "AH_5", "AH_6", "NA1", "ARC_1","ARC_2", "ARC_3", "ARC_4", "ARC_5", "ARC_6", "ARC_7", "ARC_8", "ARC_9", "ARC_10", "ARC_11", "ARC_12", "ARC_13", "ARC_14", "ARC_15", "ARC_16", "ARC_17", "ARC_18", "ARC_19", "ARC_20", "ARC_21", "ARC_22", "ARC_23", "ARC_24", "ARC_25", "ARC_26", "ARC_27", "NA2",  "DMH_1", "DMH_2", "DMH_3", "DMH_4", "NA3", "ID_1",  "ID_2", "ID_3", "ID_4", "ID_5", "ID_6", "ID_7", "ID_8", "ID_9", "ID_10", "ID_11", "ID_12", "ID_13", "ID_14","NA4", "LH_1", "LH_2", "LH_3", "LH_4", "NA5","MN_1", "MN_2", "MN_3", "MN_4", "MN_5", "MN_6", "MN_7", "MN_8", "MN_9","NA6", "PO_1", "PO_2", "PO_3", "PO_4", "PO_5", "PO_6", "PO_7", "PO_8", "NA7","PVH_1",  "PVH_2", "PVH_3", "PVH_4", "PVH_5", "PVH_6", "PVH_7", "PVH_8", "PVH_9","PVH_10", "PVH_11", "PVH_12", "PVH_13", "PVH_14", "PVH_15", "PVH_16", "PVH_17", "PVH_18", "NA8", "SCN_1", "SCN_10", "SCN_2", "SCN_3", "SCN_4", "SCN_5", "SCN_6", "SCN_7", "SCN_8", "SCN_9","NA9", "SMN_1", "SMN_10", "SMN_11", "SMN_2", "SMN_3", "SMN_4", "SMN_5", "SMN_6", "SMN_7", "SMN_8", "SMN_9","NA10",  "TM_1", "TM_2", "TM_3", "TM_4", "TM_5", "TM_6", "NA11","VMH_1","VMH_2", "VMH_3", "VMH_4", "VMH_5", "VMH_6", "VMH_7", "VMH_8", "VMH_9",  "VMH_10", "VMH_11", "VMH_12", "VMH_13", "VMH_14", "VMH_15", "VMH_16", "VMH_17", "VMH_18", "VMH_19", "VMH_20", "VMH_21", "VMH_22", "VMH_23", "VMH_24", "VMH_25",  "NA12",  "ZI_1", "ZI_2", "ZI_3", "ZI_4", "ZI_5"))

OutputHeatMap = ggplot(Hicat_gath2, aes(x = Gene, y = Cluster, fill = Log2)) +
  geom_tile()+coord_flip() + scale_fill_gradient( low = "#eeefef", high = "#2d75a4", na.value = "#FFFFFF") + theme(axis.title = element_blank(), axis.text.y = element_text(face = "italic"), axis.text.x = element_text(angle =90, hjust = 0.5, vjust = 1))

pdf("Fig2_Hicats_heatmap.pdf", width = 6, height = 4)
print(OutputHeatMap)
dev.off()

GetUMAP = as.data.frame(AdultNeuroExNac_V2@reductions$pca@cell.embeddings)
Assignw3D_Fig2 = merge(GetUMAP, Hicat_Nuclei2,by = 0)

ColPal2 = c("AH_1" = "#9E1B44","AH_2" = "#AB3B5E","AH_3" = "#B95C79","AH_4" = "#C77C94","AH_5" = "#D59DAE","AH_6" = "#E3BDC9","ARC_1" = "#F36C44", "ARC_2" = "#F3714A", "ARC_3" = "#F37651","ARC_4" = "#F47B58","ARC_5" = "#F4815E","ARC_6" = "#F58665","ARC_7" = "#F58B6C","ARC_8" = "#F59072","ARC_9" = "#F69579","ARC_10" = "#F69B80","ARC_11" = "#F7A086","ARC_12" = "#F7A58D","ARC_13" = "#F8AB94","ARC_14" = "#F8B09A","ARC_15" ="#F9B5A1" ,"ARC_16" = "#F9BAA8","ARC_17" = "#F9C0AE","ARC_18" = "#FAC5B5","ARC_19" = "#FACABC", "ARC_20" = "#FBCFC2","ARC_21" = "#FBD4C9","ARC_22" ="#FCDAD0" ,"ARC_23" ="#FCDFD6" ,"ARC_24" = "#FCE4DD","ARC_25" = "#FDE9E4","ARC_26" = "#FDEFEA","ARC_27" = "#FEF4F1","DMH_1" = "#FEE08B","DMH_2" = "#FEE6A2","DMH_3" = "#FEECB9","DMH_4" = "#FEF2D0", "ID_1" = "#ACD7A5", "ID_2" = "#B1D9AB","ID_3" = "#B7DCB1","ID_4" = "#BCDFB7","ID_5" = "#C2E1BD","ID_6" = "#C7E4C3","ID_7" = "#CDE7C9","ID_8" = "#D2E9CF","ID_9" = "#D8ECD5","ID_10" = "#DDEFDB","ID_11" = "#E3F1E1","ID_12" = "#E8F4E7","ID_13" =  "#EEF7ED" ,"ID_14" = "#F3F9F3","LH_1" = "#537A63","LH_2" = "#759482","LH_3" = "#97AFA1" ,"LH_4" = "#BAC9C0","MN_1" = "#2D75A4","MN_2" = "#4282AD","MN_3" = "#5790B6","MN_4" = "#6C9EBF","MN_5" ="#81ACC8" ,"MN_6" ="#96BAD1" ,"MN_7" ="#ABC7DA" ,"MN_8" = "#C0D5E3","MN_9" = "#D5E3EC","PO_1" = "#5EBADF","PO_2" = "#6FC1E2","PO_3" = "#81C9E6","PO_4" = "#93D1E9","PO_5" = "#A5D8ED","PO_6" = "#B7E0F0","PO_7" ="#C9E8F4" ,"PO_8" = "#DBEFF7","PVH_1" = "#81AFDD", "PVH_2" = "#87B3DE","PVH_3" = "#8EB7E0","PVH_4" = "#94BBE2","PVH_5" = "#9BBFE4","PVH_6" = "#A2C4E5","PVH_7" = "#A8C8E7","PVH_8" = "#AFCCE9","PVH_9" = "#B6D0EB","PVH_10" = "#BCD4ED","PVH_11" = "#C3D9EE","PVH_12" = "#C9DDF0","PVH_13" = "#D0E1F2","PVH_14" = "#D7E5F4","PVH_15" = "#DDE9F6","PVH_16" = "#E4EEF7","PVH_17" = "#EBF2F9","PVH_18" ="#F1F6FB" ,"SCN_1" = "#6E5AA7","SCN_2" = "#7B69AF","SCN_3" = "#8878B7","SCN_4" = "#9587BF","SCN_5" = "#A296C7","SCN_6" = "#AFA5CF","SCN_7" = "#BDB3D7","SCN_8" = "#CAC3DF","SCN_9" = "#D7D2E7","SCN_10" ="#E4E1EF" ,"SMN_1" = "#9D5BA5","SMN_2" = "#A568AC","SMN_3" = "#AD76B4","SMN_4" = "#B583BB","SMN_5" ="#BD91C3" ,"SMN_6" = "#C59FCA","SMN_7" = "#CEADD2","SMN_8" = "#D6BAD9","SMN_9" = "#DEC8E1","SMN_10" ="#E6D5E8" ,"SMN_11" = "#EEE3F0","TM_1" = "#ECBED8","TM_2" =  "#EEC7DD","TM_3" = "#F1D0E3","TM_4" = "#F4D9E8","TM_5" = "#F6E3EE","TM_6" = "#F9ECF3","VMH_1" ="#B8258F" ,"VMH_2" = "#BA2D93", "VMH_3" = "#BD3597","VMH_4" = "#C03E9B","VMH_5" = "#C246A0","VMH_6" = "#C54EA4","VMH_7" = "#C857A8","VMH_8" = "#CB5FAD","VMH_9" = "#CD68B1","VMH_10" = "#D070B5","VMH_11" = "#D378BA","VMH_12" = "#D681BE","VMH_13" = "#D889C2","VMH_14" = "#DB92C7","VMH_15" = "#DE9ACB","VMH_16" = "#E0A2CF","VMH_17" = "#E3ABD3","VMH_18" = "#E6B3D8","VMH_19" = "#E9BBDC", "VMH_20" = "#EBC4E0","VMH_21" = "#EECCE5","VMH_22" ="#F1D5E9" ,"VMH_23" = "#F4DDED","VMH_24" = "#F6E5F2","VMH_25" = "#F9EEF6","ZI_1" = "#731E5F","ZI_2" = "#8A4379","ZI_3" = "#A16994","ZI_4" = "#B98EAF","ZI_5" = "#D0B3C9")



set.seed(45)
StartTime = "CS22"
set.kparam =  c(40)
kval = 40
SeuName = "AdultNeuroExNac_seed45"
get.PCs = 20
get.UMAPs = 10
get.spread = 10
get.node = 1
Filename = paste(SeuName, get.PCs, get.UMAPs, get.spread, sep="_")
AdultNeuroExNac_V2@meta.data$sampletype = gsub("_.*", "",AdultNeuroExNac_V2@meta.data$sample)
AdultNeuroExNac_V2@meta.data$sampletype = gsub("HTHpo", "Adult_PO",  gsub("HTHtub", "Adult_Tub", gsub("HTHmn", "Adult_MN", gsub("HTHso", "Adult_SO",gsub("22T", "22", AdultNeuroExNac_V2@meta.data$sampletype)))))

AdultNeuroExNac_V2@meta.data$Timepoint = gsub("_.*", "", gsub("HTHpo", "Adult",  gsub("HTHtub", "Adult", gsub("HTHmn", "Adult", gsub("HTHso", "Adult",  AdultNeuroExNac_V2@meta.data$sample)))))
DefaultAssay(AdultNeuroExNac_V2) = "integrated"
AdultNeuroExNac_V2 = RunPCA(AdultNeuroExNac_V2, npcs = get.PCs)
AdultNeuroExNac_V2 <- RunUMAP(AdultNeuroExNac_V2, dims = 1:get.UMAPs, spread= get.spread, n.components = 3)
DefaultAssay(AdultNeuroExNac_V2) = "RNA"
### monocle 3
geneMeta = data.frame(gene_short_name=rownames(AdultNeuroExNac_V2))    
rownames(geneMeta) = rownames(AdultNeuroExNac_V2)
DefaultAssay(AdultNeuroExNac_V2) = 'integrated' 
redSeu = subset(AdultNeuroExNac_V2,features=rownames(geneMeta))
dat = as.matrix(redSeu@assays$RNA@counts)
cmeta = AdultNeuroExNac_V2@meta.data
#cmeta$Hicat_Nuclei2 = NULL
M3Seu <- new_cell_data_set(dat, cell_metadata = cmeta, gene_metadata = geneMeta)
M3Seu <- preprocess_cds(M3Seu, num_dim = get.PCs)
M3Seu <- reduce_dimension(M3Seu)
M3Seu@preprocess_aux$gene_loadings = redSeu@reductions[["pca"]]@feature.loadings[,1:get.PCs] ## PC by gene - replace Monocle pc's with seurat calculated
reducedDims(M3Seu)[['UMAP']] = redSeu@reductions[["umap"]]@cell.embeddings ## from code - UMAP by cell 
Filename = paste(SeuName, get.PCs, get.UMAPs, get.spread, kval, sep="_")
dir.create(file.path('~/Dropbox/Columbia/Fetal scRNAseq/Science Paper - 2022/Neuronal_Optimization/', Filename))

M3Seu <- cluster_cells(M3Seu,reduction_method = "UMAP", k = kval)
M3Seu <- learn_graph(M3Seu)

Pull = as.data.frame(paste('Y_',M3Seu@principal_graph_aux[['UMAP']]$pr_graph_cell_proj_closest_vertex[,1],sep=''))
row.names(Pull) = colnames(M3Seu)
colnames(Pull) = "Vertex"
Pull$Timepoint = gsub("_.*", "", gsub("HTHpo", "Adult",  gsub("HTHtub", "Adult", gsub("HTHmn", "Adult", gsub("HTHso", "Adult", row.names(Pull))))))
Pull$Val = 1
Pull_n = as.data.frame(table(Pull$Timepoint))
Pull3 = Pull %>% dplyr::group_by(Vertex, Timepoint) %>% dplyr::summarise(n_ = sum(Val, na.rm=T))
Pull4 = as.data.frame(table(Pull$Vertex))
Pull5 = merge(Pull3, Pull4, by.x = "Vertex", by.y = "Var1")
Pull5$Percent = Pull5$n_/Pull5$Freq
PullCS22 = subset(Pull5, Pull5$Timepoint == StartTime & Pull5$Freq > 10)
PullCS22 = PullCS22[order(-PullCS22$Percent), ]

Pull6 = subset(Pull, Pull$Vertex %in% PullCS22[[1]][get.node])

M3SeuNODE = M3Seu
M3SeuNODE <- order_cells(M3SeuNODE, root_pr_nodes=PullCS22[[1]][get.node])
dir.create(file.path(paste('~/Dropbox/Columbia/Fetal scRNAseq/Science Paper - 2022/Neuronal_Optimization/', Filename, sep=""), paste("Node", get.node, "_", PullCS22[[1]][get.node], sep="")))

############## FIGURE 2B ############## 
p = plot_cells_3d(M3SeuNODE,color_cells_by = "Timepoint", color_palette = c("#3288bd", "#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#e6f598", "#abdda4", "#66c2a5"))
htmlwidgets::saveWidget(p, paste('~/Dropbox/Columbia/Fetal scRNAseq/Science Paper - 2022/Neuronal_Optimization/', Filename, "/Node", get.node, "_", PullCS22[[1]][get.node], '/', Filename, "K_", kval, "_Sample.html", sep=""))

############## FIGURE 2C ############## 
p = plot_cells_3d(M3SeuNODE,color_cells_by = "Hicat_Nuclei2", color_palette = ColPal2)
htmlwidgets::saveWidget(p, paste('~/Dropbox/Columbia/Fetal scRNAseq/Science Paper - 2022/Neuronal_Optimization/', Filename, "/Node", get.node, "_", PullCS22[[1]][get.node], '/', Filename, "K_", kval, "_Hicat_Nuclei2.html", sep=""))

############## FIGURE 2E ############## 
p = plot_cells_3d(M3SeuNODE,color_cells_by = "pseudotime")
htmlwidgets::saveWidget(p, paste('~/Dropbox/Columbia/Fetal scRNAseq/Science Paper - 2022/Neuronal_Optimization/', Filename, "/Node", get.node, "_", PullCS22[[1]][get.node], '/', Filename, "K_", kval, "_Pseudotime.html", sep=""))


p = plot_cells_3d(M3SeuNODE,color_cells_by = "partition")
htmlwidgets::saveWidget(p, paste('~/Dropbox/Columbia/Fetal scRNAseq/Science Paper - 2022/Neuronal_Optimization/', Filename, "/Node", get.node, "_", PullCS22[[1]][get.node], '/', Filename,"K_", kval, "_Partition.html", sep="")) #Want to check that the majority of the cells are contained in one partition
tmpLin = M3SeuNODE@principal_graph_aux[['UMAP']]["pseudotime"][[1]]
tmpLin[which(tmpLin==Inf)] =50

############## FIGURE 2F ############## 
M3SeuNODEpseu = tapply(X=tmpLin,INDEX=Pull$Vertex,FUN=mean, na.rm = TRUE)
mleaves = monocle3:::leaf_nodes(M3Seu,reduction_method='UMAP') 
colData(M3Seu)$MonocleLeaf = NA
colData(M3Seu)$MonocleLeaf[colData(M3Seu)$Vertex%in%names(mleaves)] = colData(M3Seu)$Vertex[colData(M3Seu)$Vertex%in%names(mleaves)]
p=plot_cells_3d(M3Seu,color_cells_by = "MonocleLeaf",color_palette=hue_pal()(length(unique(colData(M3Seu)$MonocleLeaf))))
htmlwidgets::saveWidget(p, paste('~/Dropbox/Columbia/Fetal scRNAseq/Science Paper - 2022/Neuronal_Optimization/', Filename, "/Node", get.node, "_", PullCS22[[1]][get.node], '/', Filename, "K_", kval, "_MonocleLeaves.html", sep=""))


mbranches = monocle3:::branch_nodes(M3Seu,reduction_method='UMAP') 
colData(M3Seu)$MonocleBranch = NA
colData(M3Seu)$MonocleBranch[colData(M3Seu)$Vertex%in%names(mbranches)] = colData(M3Seu)$Vertex[colData(M3Seu)$Vertex%in%names(mbranches)]
p=plot_cells_3d(M3Seu,color_cells_by = "MonocleBranch",color_palette=hue_pal()(length(unique(colData(M3Seu)$MonocleBranch))))
htmlwidgets::saveWidget(p, paste('~/Dropbox/Columbia/Fetal scRNAseq/Science Paper - 2022/Neuronal_Optimization/', Filename, "/Node", get.node, "_", PullCS22[[1]][get.node], '/', Filename, "K_", kval, "_MonocleBranches.html", sep=""))


colData(M3Seu)$Leaves_Branches = NA
colData(M3Seu)$Leaves_Branches[colData(M3Seu)$Vertex%in%c(names(mleaves), names(mbranches))] = colData(M3Seu)$Vertex[colData(M3Seu)$Vertex%in%c(names(mleaves), names(mbranches))]

colrange_branches <- colorRampPalette(c("#d1e3d1", "#087508")) #green
colrange_branches_n <- colrange_branches(length(unique(mbranches)))
colrange_branches_ndf = as.data.frame(colrange_branches_n)
colrange_branches_ndf$Vertex = unique(mbranches)
colnames(colrange_branches_ndf) = c("Colour", "Vertex")

colrange_leaves <- colorRampPalette(c("#ffd4b3", "#f77007")) #orange
colrange_leaves_n <- colrange_leaves(length(unique(mleaves)))
colrange_leaves_ndf = as.data.frame(colrange_leaves_n)
colrange_leaves_ndf$Vertex = unique(mleaves)
colnames(colrange_leaves_ndf) = c("Colour", "Vertex")

CombineColors = rbind(colrange_leaves_ndf, colrange_branches_ndf)
CombineColors = CombineColors[order(CombineColors$Vertex), ]
OutsColor =  CombineColors$Colour
names(OutsColor) = paste("Y_", CombineColors$Vertex, sep="")

p=plot_cells_3d(M3Seu,color_cells_by = "Leaves_Branches",color_palette=OutsColor)
htmlwidgets::saveWidget(p, paste('~/Dropbox/Columbia/Fetal scRNAseq/Science Paper - 2022/Neuronal_Optimization/', Filename, "/Node", get.node, "_", PullCS22[[1]][get.node], '/', Filename, "K_", kval, "_MonocleLeavesandBranches.html", sep=""))  #### FOR FIG 2 iv #### 


############## FIGURE 2G ############## 
##Generate Lineage Plots
rtti287m = data.frame(root=rep(PullCS22[[1]][get.node],length(unique(mleaves))), leaf=paste("Y_", unique(mleaves), sep=""))
PullPartitions = as.data.frame(partitions(M3Seu, reduction_method = "UMAP"))
colnames(PullPartitions) = "Partitions"
TablePartitions = as.data.frame(table(PullPartitions$Partitions))
VertexInPartit = subset(PullPartitions, row.names(PullPartitions) %in% row.names(Pull6))
OtherVertexinPart = subset(PullPartitions, PullPartitions$Partitions %in% unique(VertexInPartit$Partitions))
OtherVertexinPull = subset(Pull, row.names(Pull) %in% row.names(OtherVertexinPart))

rtti287m2 = subset(rtti287m, rtti287m$leaf %in% OtherVertexinPull$Vertex)


linksNODE = data.frame(matrix(ncol = 3, nrow = 0))
colnames(linksNODE) = c("from", "to", "weight")
for(i in seq(1, nrow(rtti287m2), 1)){
FROM=rtti287m2$root[i]
TO=rtti287m2$leaf[i]
lins = all_simple_paths(M3Seu@principal_graph[['UMAP']],FROM,TO)
lengs = unlist(lapply(lins,FUN=length))
shortInd = which(lengs==min(lengs))
if(length(lins) > 0){
tmpLin = lins[[shortInd]]
tmpEdge = data.frame(from = paste('Y_',as.vector(tmpLin[-length(tmpLin)]),sep=''),to=paste('Y_',as.vector(tmpLin[-1]),sep=''))
tmpWeight = abs(M3SeuNODEpseu[tmpEdge$to] - M3SeuNODEpseu[tmpEdge$from])
tmpEdge$weight = tmpWeight
if(dim(linksNODE)[1]==0){
linksNODE = tmpEdge
} else {
linksNODE = rbind(linksNODE,tmpEdge)
}
}
}

linksNODE = unique(linksNODE)
nodesNODE = data.frame(id=unique(c(linksNODE$from,linksNODE$to)))
nodesNODE$type = 'lineage'
nodesNODE$type[nodesNODE$id%in%names(table(linksNODE$from)[which(table(linksNODE$from)>1)])] = 'branch' 
nodesNODE$type[nodesNODE$id==PullCS22[[1]][get.node]] = 'root'
nodesNODE$type[nodesNODE$id%in%names(table(c(linksNODE$from,linksNODE$to))[which(table(c(linksNODE$from,linksNODE$to))==1)])] = 'leaf'
nodesNODE$type2 = 4
nodesNODE$type2[nodesNODE$id%in%names(table(linksNODE$from)[which(table(linksNODE$from)>1)])] = 3
nodesNODE$type2[nodesNODE$id==PullCS22[[1]][get.node]] = 1
nodesNODE$type2[nodesNODE$id%in%names(table(c(linksNODE$from,linksNODE$to))[which(table(c(linksNODE$from,linksNODE$to))==1)])] = 2
nodesNODE$rbl = nodesNODE$id
nodesNODE$rbl[nodesNODE$type%in%'lineage']='' ## can add other columns with metadata
net.NODE <- graph_from_data_frame(d=linksNODE, vertices=nodesNODE, directed=T)
V(net.NODE)$size <- c(10,5,3,2)[V(net.NODE)$type2]
V(net.NODE)$frame.color <- "white"
V(net.NODE)$color <- c("darkgray","skyblue3","skyblue2","gray")[V(net.NODE)$type2]
V(net.NODE)$label <- V(net.NODE)$rbl #V(net.NODE)$label <- ""
E(net.NODE)$arrow.mode <- 0
lrt.NODE = layout.reingold.tilford(net.NODE)
#End Generation of Lineage Plots

pdf(paste('~/Dropbox/Columbia/Fetal scRNAseq/Science Paper - 2022/Neuronal_Optimization/', Filename, "/Node", get.node, "_", PullCS22[[1]][get.node], '/', Filename, "K_", kval, "ShortestPath_MainLineage.pdf", sep=""),width=12,height=8)
print(plot(net.NODE,layout=lrt.NODE,vertex.label.cex=0.5))
dev.off() #Save version with labels at branches and leaved

V(net.NODE)$label <- ""
pdf(paste('~/Dropbox/Columbia/Fetal scRNAseq/Science Paper - 2022/Neuronal_Optimization/', Filename, "/Node", get.node, "_", PullCS22[[1]][get.node], '/', Filename, "K_", kval, "ShortestPath_MainLineageNoLabs.pdf", sep=""),width=12,height=8)
print(plot(net.NODE,layout=lrt.NODE,vertex.label.cex=0.5))
dev.off() #Save version without labels - used for Fig 2vi

V(net.NODE)$label <- V(net.NODE)$name #V(net.NODE)$label <- ""
pdf(paste('~/Dropbox/Columbia/Fetal scRNAseq/Science Paper - 2022/Neuronal_Optimization/', Filename, "/Node", get.node, "_", PullCS22[[1]][get.node], '/', Filename, "K_", kval, "ShortestPath_MainLineageWithLabels.pdf", sep=""),width=12,height=8)
print(plot(net.NODE,layout=lrt.NODE,vertex.label.cex=0.5))
dev.off() #Save version with ALL labels - used to assign nuclei


SeuFile2 = AdultNeuroExNac_V2
PullVertex = Pull %>% select(Vertex)
SeuFile2= AddMetaData(SeuFile2, PullVertex, "Vertex")
ind = which(SeuFile2@meta.data$Vertex%in%nodesNODE$id)

############## FIGURE 2H + SUPP ############## 
#Generate Lineages with gene expression overlaid - included in Figure 2v and in supplementals
for(i in c("DLX1", "DLX2", "IRX3", "KISS1")){
gene = i
tmpGene = SeuFile2@assays$RNA@counts[gene,ind]
tmpGeneVertex = tapply(X=tmpGene,INDEX=SeuFile2@meta.data$Vertex[ind],FUN=mean)
tmpGeneVertex2 = round(scales::rescale(tmpGeneVertex, to = c(1, 100)),0)
fun_color_range <- colorRampPalette(c("#919191", "red"))
my_colors <- fun_color_range(100)
tmpGeneColors = my_colors[tmpGeneVertex2]
names(tmpGeneColors) = names(tmpGeneVertex2)
V(net.NODE)$color <- tmpGeneColors[as.vector(nodesNODE$id)]

pdf(paste('~/Dropbox/Columbia/Fetal scRNAseq/Science Paper - 2022/Neuronal_Optimization/', Filename, "/Node", get.node, "_", PullCS22[[1]][get.node], '/', Filename, "K_", kval, "_ShortestPath_Lineage_", i, ".pdf", sep=""),width=12,height=8)
print(plot(net.NODE,layout=lrt.NODE, vertex.label =NA, label = i))
dev.off()
}


GenerateMetaDataVertex = function(ListMeta){
NodesforDEGs = as.data.frame(matrix(ncol =1, nrow=0))
colnames(NodesforDEGs) = c("Node")
for(x in names(SplitNodes)){
GetNodes = SplitNodes[[x]]
Pull_Spec = subset(Pull, Pull$Vertex %in% GetNodes)  
Pull_Spec$Node = x
Pull_Spec2 = Pull_Spec %>% dplyr::select("Node")
NodesforDEGs = rbind(NodesforDEGs, Pull_Spec2)
}
return(NodesforDEGs)
}

#### EMBRYO + ADULT #### NUCLEI ASSIGNS
#load("~/FETAL_FIG2_ADULTEMBRYODATA.RData")
SplitNodes = list(Intermediates = c("Y_8", "Y_74", "Y_137", "Y_164", "Y_39", "Y_48", "Y_163",  "Y_68,", "Y_85", "Y_21", "Y_78", "Y_203", "Y_38", "Y_241", "Y_171", "Y_81", "Y_172", "Y_26", "Y_66", "Y_258", "Y_40", "Y_73", "Y_95", "Y_100", "Y_101", "Y_65", "Y_25", "Y_75", "Y_61", "Y_62", "Y_27", "Y_120", "Y_97", "Y_80", "Y_46", "Y_365", "Y_237", "Y_377", "Y_47", "Y_33", "Y_44", "Y_82", "Y_29", "Y_110", "Y_130", "Y_364", "Y_363"), ARC = c( "Y_98", "Y_77", "Y_79", "Y_90", "Y_58", "Y_202", "Y_63", "Y_132", "Y_51", "Y_382", "Y_379", "Y_261", "Y_189", "Y_10", "Y_213"), VMH = c("Y_168", "Y_225", "Y_249", "Y_234", "Y_251", "Y_239", "Y_229", "Y_232", "Y_373", "Y_360", "Y_238", "Y_226", "Y_242", "Y_217", "Y_247", "Y_106", "Y_105", "Y_157", "Y_187", "Y_375", "Y_208", "Y_153", "Y_184", "Y_60", "Y_59", "Y_369", "Y_16", "Y_111", "Y_64", "Y_72", "Y_165", "Y_175", "Y_83", "Y_5", "Y_52", "Y_31", "Y_57", "Y_107", "Y_70", "Y_87", "Y_374", "Y_158", "Y389"), PVH = c("Y_245", "Y_99", "Y_1", "Y_291", "Y_3", "Y_7", "Y_109", "Y_12", "Y_311", "Y_20", "Y_9", "Y_19", "Y_88", "Y_129", "Y_34", "Y_93", "Y_122", "Y_28", "Y_174", "Y_169"), SMN = c("Y_322", "Y_350", "Y_115", "Y_359", "Y_145", "Y_181", "Y_339", "Y_162", "Y_178", "Y_269", "Y_142", "Y_149", "Y_67", "Y_54", "Y_126", "Y_240", "Y_32", "Y_214", "Y_182", "Y_362", "Y_53", "Y_315", "Y_192", "Y_356", "Y_18", "Y_112", "Y_215", "Y_4", "Y_387", "Y_45", "Y_199"), TM = c("Y_37", "Y_15", "Y_69", "Y_150", "Y_55", "Y_91", "Y_102", "Y_108", "Y_167", "Y_190", "Y_84", "Y_56", "Y_116"), SCN = c("Y_49", "Y_86", "Y_119", "Y_17", "Y_104", "Y_207", "Y_118", "Y_136", "Y_180", "Y_152", "Y_156", "Y_159", "Y_196", "Y_151", "Y_385", "Y_250", "Y_179", "Y_147"),  MN = c( "Y_296", "Y_278", "Y_354", "Y_260", "Y_346", "Y_334", "Y_277", "Y_263", "Y_351", "Y_206", "Y_274", "Y_305", "Y_320", "Y_280", "Y_357", "Y_299", "Y_275", "Y_155", "Y_319", "Y_177", "Y_330", "Y_185", "Y_332", "Y_314", "Y_337", "Y_188", "Y_261", "Y_293", "Y_212", "Y_340", "Y_279", "Y_348", "Y_176", "Y_276", "Y_325", "Y_270", "Y_331", "Y_160", "Y_266", "Y_391", "Y_378", "Y_220", "Y_273", "Y_345", "Y_204", "Y_198", "Y_117", "Y_113", "Y_281", "Y_23", "Y_223", "Y_231", "Y_139", "Y_35", "Y_166", "Y_134", "Y_141", "Y_22", "Y_290", "Y_92", "Y_94", "Y_140", "Y_133", "Y_367", "Y_96", "Y_347", "Y_259", "Y_333", "Y_201", "Y_318", "Y_161", "Y_298", "Y_284", "Y_282", "Y_295", "Y_265", "Y_310", "Y_285", "Y_327", "Y_309", "Y_342", "Y_343", "Y_211", "Y_349", "Y_308", "Y_170", "Y_326", "Y_264", "Y_335", "Y_262"), ID = c("Y_135", "Y_252", "Y_218", "Y_254", "Y_236", "Y_390", "Y_244", "Y_256", "Y_224", "Y_216", "Y_233", "Y_186", "Y_127", "Y_248", "Y_131", "Y_228", "Y_230", "Y_221", "Y_121", "Y_30", "Y_41", "Y_191", "Y_50", "Y_103", "Y_324", "Y_257", "Y_255", "Y_13", "Y_227", "Y_222", "Y_253", "Y_386", "Y_243", "Y_246"), LH = c("Y_42", "Y_71", "Y_114", "Y_123", "Y_11", "Y_384", "Y_205", "Y_371", "Y_6", "Y_24", "Y_173", "Y_89", "Y_235", "Y_193", "Y_76"), ZI = c("Y_353", "Y_301", "Y_336", "Y_316", "Y_344", "Y_303", "Y_267", "Y_317", "Y_341", "Y_355", "Y_283"))



EmbryoAdultNuclei = GenerateMetaDataVertex(SplitNodes)
AdultNeuroExNac_V3 = subset(AdultNeuroExNac_V2, cells = row.names(EmbryoAdultNuclei))

AdultNeuroExNac_V3 = AddMetaData(AdultNeuroExNac_V3, EmbryoAdultNuclei, "EmbryoAdultNuclei")
Idents(AdultNeuroExNac_V3) = "EmbryoAdultNuclei"
DimPlot(AdultNeuroExNac_V3, label =T) +NoLegend()

AdultNeuroExNac_V3@meta.data$Timepoint = gsub("_.*", "",gsub("GW22T", "GW22", gsub("HTHpo", "Adult",  gsub("HTHtub", "Adult", gsub("HTHmn", "Adult", gsub("HTHso", "Adult", AdultNeuroExNac_V3@meta.data$sample))))))

AdultNeuroExNac_V3@meta.data$Trimester = ifelse(AdultNeuroExNac_V3@meta.data$Timepoint %in% c("CS22"), "Tri1", AdultNeuroExNac_V3@meta.data$Timepoint)
AdultNeuroExNac_V3@meta.data$Trimester = ifelse(AdultNeuroExNac_V3@meta.data$Trimester %in% c("GW16", "GW18",  "GW19", "GW20",  "GW22", "GW25"), "Tri2", AdultNeuroExNac_V3@meta.data$Trimester)
AdultNeuroExNac_V3@meta.data$TriXNode = paste(AdultNeuroExNac_V3@meta.data$EmbryoAdultNuclei, AdultNeuroExNac_V3@meta.data$Trimester, sep=" - ")

AdultNeuroExNac_V3@meta.data$TriXNode = factor(AdultNeuroExNac_V3@meta.data$TriXNode, levels = c("TM - Tri1", "TM - Tri2","TM - Adult", "ARC - Tri1", "ARC - Tri2", "ARC - Adult", "PVH - Tri1", "PVH - Tri2", "PVH - Adult",   "VMH - Tri1", "VMH - Tri2", "VMH - Adult", "LH - Tri1" , "LH - Tri2", "LH - Adult" , "DMH - Tri1", "DMH - Tri2", "DMH - Adult","SCN - Tri1", "SCN - Tri2","SCN - Adult", "SMN - Tri1", "SMN - Tri2", "SMN - Adult", "MN - Tri1", "MN - Tri2", "MN - Adult", "ZI - Tri1", "ZI - Tri2", "ZI - Adult", "PO - Tri1", "PO - Tri2","PO - Adult","ID - Tri1","ID - Tri2","ID - Adult", "Intermediates - Tri1", "Intermediates - Tri2", "Intermediates - Adult"))

sz=8

############## FIGURE 2I ############## 
Idents(AdultNeuroExNac_V3) = "TriXNode"
DefaultAssay(AdultNeuroExNac_V3) = "RNA"
DoPl = DotPlot(AdultNeuroExNac_V3, features = unique(rev(c("HDC", "TBX3", "POMC","HMX2","GHRH", "GAL",  "AGRP", "NPY",   "SIM1", "POU3F2","OXT", "AVP", "CRH", 'NPTX2', "FEZF1", "NR5A1", "SOX14",  "LHX9", "LHX8", "LHX9", "HCRT","SIX3", "SIX6", "BARHL1", "LMX1A", "IRX3", "IRX5", "FOXA1", "FOXB1", "LHX1", "PITX2", "ARX", "LHX6"))), assay="RNA", cols=c("lightgrey", "navy"), dot.min = 0.02)+coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size=sz*1), axis.text.y = element_text(size=sz*1, face = "italic"), legend.text = element_text(size=sz), legend.title = element_text(size=sz)) +xlab("")+ylab("")

pdf(paste("Fetal_KeyGenes_DotPlot_EmbryoAdult.pdf", sep=""), width = 6, height =6) 
print(DoPl)
dev.off() 
#, dot.scale = 4.5


#### EMBRYO + ADULT #### TEMPORAL ASSIGNS
############## FIGURE 2J ############## 
SplitNodes = list(Node1 = c("Y_8"), Node2 = c("Y_74", "Y_137", "Y_164", "Y_39", "Y_48", "Y_163"), Node3 = c("Y_37", "Y_15", "Y_69", "Y_150", "Y_55", "Y_91", "Y_102", "Y_108", "Y_167", "Y_190", "Y_84", "Y_56", "Y_116"), Node4 = c("Y_68,", "Y_85", "Y_21", "Y_78", "Y_203", "Y_38", "Y_241", "Y_171", "Y_81"), Node5 = c("Y_98", "Y_77", "Y_79", "Y_90", "Y_58", "Y_202", "Y_63"), Node6 = c("Y_172", "Y_26", "Y_66", "Y_258"), Node7 = c("Y_40", "Y_73", "Y_95", "Y_100", "Y_101", "Y_65", "Y_25", "Y_75", "Y_61", "Y_62", "Y_27", "Y_120", "Y_97", "Y_80", "Y_46", "Y_121", "Y_30", "Y_41", "Y_191", "Y_50", "Y_103", "Y_324", "Y_257", "Y_255", "Y_13", "Y_227", "Y_222", "Y_253", "Y_386", "Y_243", "Y_246"), Node8 = c("Y_149", "Y_67", "Y_54", "Y_126", "Y_240", "Y_32", "Y_214", "Y_182", "Y_362", "Y_53", "Y_315", "Y_192", "Y_356", "Y_18", "Y_112", "Y_215", "Y_4", "Y_387", "Y_45", "Y_199"), Node9 = c("Y_245", "Y_99", "Y_1", "Y_291", "Y_3", "Y_7", "Y_109", "Y_12", "Y_311", "Y_20", "Y_9", "Y_19", "Y_88", "Y_129", "Y_34", "Y_93", "Y_122", "Y_28", "Y_174", "Y_169"), Node10 = c("Y_49", "Y_86", "Y_119", "Y_17", "Y_104", "Y_207", "Y_118", "Y_136", "Y_180", "Y_152", "Y_156", "Y_159", "Y_196", "Y_151", "Y_385", "Y_250", "Y_179", "Y_147"),  Node11 = c("Y_365", "Y_237", "Y_377", "Y_47", "Y_33", "Y_44", "Y_82", "Y_29", "Y_110"), Node12 = c("Y_322", "Y_350", "Y_115", "Y_359", "Y_145", "Y_181", "Y_339", "Y_162", "Y_178", "Y_269", "Y_142"), Node13 = c("Y_130", "Y_364", "Y_363", "Y_132", "Y_51", "Y_382", "Y_379", "Y_261", "Y_189", "Y_10", "Y_213"), Node14 = c("Y_42", "Y_71", "Y_114", "Y_123", "Y_11", "Y_384", "Y_205", "Y_371", "Y_6", "Y_24", "Y_173", "Y_89", "Y_235", "Y_193", "Y_76"), Node15 = c("Y_31", "Y_57", "Y_107", "Y_70", "Y_87", "Y_374", "Y_158", "Y389"), Node16 = c("Y_106", "Y_105", "Y_157", "Y_187", "Y_375", "Y_208", "Y_153", "Y_184", "Y_60", "Y_59", "Y_369", "Y_16", "Y_111", "Y_64", "Y_72", "Y_165", "Y_175", "Y_83", "Y_5", "Y_52"), Node17 = c("Y_168", "Y_225", "Y_249", "Y_234", "Y_251", "Y_239", "Y_229", "Y_232", "Y_373", "Y_360", "Y_238", "Y_226", "Y_242", "Y_217", "Y_247"), Node18 = c("Y_198", "Y_117", "Y_113", "Y_281", "Y_23", "Y_223", "Y_231", "Y_139", "Y_35", "Y_166", "Y_134", "Y_141", "Y_22", "Y_290", "Y_92", "Y_94", "Y_140", "Y_133", "Y_367", "Y_96", "Y_347", "Y_259", "Y_333"), Node19 = c("Y_296", "Y_278", "Y_354", "Y_260", "Y_346", "Y_334", "Y_277", "Y_263", "Y_353", "Y_301", "Y_336", "Y_316", "Y_344", "Y_303", "Y_267", "Y_317", "Y_341", "Y_355", "Y_283", "Y_201", "Y_318", "Y_161", "Y_298", "Y_284"), Node20 = c("Y_282", "Y_295", "Y_265", "Y_310", "Y_285", "Y_327", "Y_309", "Y_342", "Y_343", "Y_211", "Y_349", "Y_308", "Y_170", "Y_326", "Y_264", "Y_335", "Y_262"), Node21 = c("Y_275", "Y_155", "Y_319", "Y_177", "Y_330", "Y_185", "Y_332", "Y_351", "Y_206", "Y_274", "Y_305", "Y_320", "Y_280", "Y_357", "Y_299"), Node22= c("Y_135", "Y_252", "Y_218", "Y_254", "Y_236", "Y_390", "Y_244", "Y_256", "Y_224", "Y_216"), Node23 = c("Y_233", "Y_186", "Y_127", "Y_248", "Y_131", "Y_228", "Y_230", "Y_221"), Node24 = c("Y_314", "Y_337", "Y_188", "Y_261", "Y_293", "Y_212", "Y_340", "Y_279", "Y_348", "Y_176", "Y_276", "Y_325", "Y_270", "Y_331", "Y_160", "Y_266", "Y_391", "Y_378", "Y_220", "Y_273", "Y_345", "Y_204"))
  

Pull$Time2 = gsub("_.*", "", row.names(Pull))
Pull$Time2 = gsub("HTHmn", "Adult_MN", gsub("HTHpo", "Adult_PO", gsub("HTHso",  "Adult_SO", gsub( "HTHtub", "Adult_TUB", Pull$Time2))))

N_Timepoint = as.data.frame(table(Pull$Time2))
CompilePercents = as.data.frame(matrix(ncol =6, nrow=0))
colnames(CompilePercents) = c("Timepoint", "N_Node", "N_All", "Percent_Time", "Percent_TimeNode", "Node")

for(x in names(SplitNodes)){
GetNodes = SplitNodes[[x]]
Pull_Spec = subset(Pull, Pull$Vertex %in% GetNodes)  
ReduceSpec = as.data.frame(table(Pull_Spec$Time2))
ReduceSpec_merge = merge(ReduceSpec, N_Timepoint, by = "Var1", all = T)
ReduceSpec_merge[is.na(ReduceSpec_merge)] <- 0
ReduceSpec_merge$Percent = ReduceSpec_merge$Freq.x / ReduceSpec_merge$Freq.y
ReduceSpec_merge$Percent2 = ReduceSpec_merge$Percent/sum(ReduceSpec_merge$Percent)*100
colnames(ReduceSpec_merge) = c("Timepoint", "N_Node", "N_All", "Percent_Time", "Percent_TimeNode")
ReduceSpec_merge$Node = x
CompilePercents = rbind(CompilePercents, ReduceSpec_merge)
}
CompilePercents$Node = gsub("Node", "", CompilePercents$Node)
CompilePercents$Node = factor(CompilePercents$Node, levels = c("1", "2", "3",  "4",  "5",  "6",  "7",  "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24"))


CompilePercents$Timepoint = gsub("22T", "22", CompilePercents$Timepoint)
CompilePercents$Trimester = ifelse(CompilePercents$Timepoint %in% c("CS13", "CS14", "CS15", "CS22"), "Tri1", as.character(CompilePercents$Timepoint))
CompilePercents$Trimester = ifelse(CompilePercents$Trimester %in% c("GW16", "GW18", "GW19", "GW20", "GW22", "GW25"), "Tri2", as.character(CompilePercents$Trimester))
CompilePercents$Trimester = ifelse(CompilePercents$Trimester %in% c("Adult_PO", "Adult_SO", "Adult_TUB", "Adult_MN"), "Adult", as.character(CompilePercents$Trimester))

CompilePercents$Timepoint = factor(CompilePercents$Timepoint, levels = c("CS13", "CS14", "CS15", "CS22", "GW16", "GW18", "GW19", "GW20", "GW22", "GW25", "Adult_PO", "Adult_SO", "Adult_TUB", "Adult_MN"))
CompilePercents$Trimester = factor(CompilePercents$Trimester, levels = c("Tri1", "Tri2", "Adult"))


############## FIGURE 2K ############## 
SP = ggplot(CompilePercents, aes(fill=Timepoint, y=Percent_TimeNode, x=Node))+ theme_classic() +
      geom_bar(position="stack", stat="identity") + theme(legend.position = "bottom", axis.text = element_text(size = 6), axis.title = element_text(size = 8)) + ylab("% cells") + xlab("")+ scale_y_continuous(expand= c(0,0))

pdf("FetalEMBRYOADULT_Fig2_Temporal_StackedbyPercent_TimepointDual_6JUN22.pdf", width = 3.7, height = 2.5)
print(SP)
dev.off()


############## FIGURE 2L ############## 
SP = ggplot(CompilePercents, aes(fill=Trimester, y=Percent_TimeNode, x=Node))+ theme_classic() +
      geom_bar(position="stack", stat="identity") + theme(legend.position = "bottom", axis.text = element_text(size = 6), axis.title = element_text(size = 8)) + ylab("% cells") + xlab("")+ scale_y_continuous(expand= c(0,0))

pdf("FetalEMBRYOADULT_Fig2_Temporal_StackedbyPercent_TrimesterDual_6JUN22.pdf", width = 3.7, height = 2)
print(SP)
dev.off()


###

COMPARE = merge(EmbryoAdultNuclei, Hicats_Nuclei, by = 0)
COMPARE$Hicat_NucleBroad = gsub("_.*", "", COMPARE$Hicat_Nuclei)

COMPARE_V2 = subset(COMPARE, ! COMPARE$Hicat_NucleBroad %in% c("Excl (< 10 cells)", "Ukn") & ! COMPARE$Node %in% "Intermediates")
COMPARE_V2$Match = COMPARE_V2$Hicat_NucleBroad == COMPARE_V2$Node#
COMPARE_TABLE = table(COMPARE_V2$Match)
COMPARE_TABLE[["TRUE"]]/(COMPARE_TABLE[["TRUE"]]+COMPARE_TABLE[["FALSE"]])

#0.7463897

#### SUPP Tables
AdultNeuroExNac_V2 = AddMetaData(AdultNeuroExNac_V2, EmbryoAdultNuclei, "EmbryoAdultNuclei")
AdultNeuroExNac_V2@meta.data$count = 1
SupTable7 = AdultNeuroExNac_V2@meta.data %>% group_by(EmbryoAdultNuclei, sampletype) %>% dplyr::summarise("nCells" = sum(count),	"nCount_mean" = mean(nCount_RNA),	"nCount_median" = median(nCount_RNA),	"nCount_sd" = sd(nCount_RNA),	"nFeature_mean" = mean(nFeature_RNA),	"nFeature_median" = median(nFeature_RNA),	"nFeature_sd" = sd(nFeature_RNA),	"percMT_mean" = mean(percent.mt),	"percMT_median"= median(percent.mt),	"percMT_sd"= sd(percent.mt),)
write.csv(SupTable7, "SupTable7.csv")

write.csv(EmbryoAdultNuclei, "SupTable8.csv")


## EmbryoAdultNuclei object saved here: 'https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/Fig2_Neurons_FINAL.RData'
