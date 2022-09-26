## to replicate conda environment used in paper, install conda, and then call the following in terminal: conda env create -f environment_droplet.yml - yml file here: https://github.com/brianherb/HumanHypothalamusDev/blob/main/environment_droplet.yml

library(dendextend)
library(matrixStats)
library(Matrix)
library(scrattch.hicat)
library(dplyr)
library(Seurat)
library(dplyr)
library(ggplot2)
library(scales)

## following tutorial at https://taxonomy.shinyapps.io/scrattch_tutorial

devtools::source_url('https://github.com/brianherb/HumanHypothalamusDev/blob/main/HypoPub_functions_reference.R?raw=TRUE')

## load in seurat object 

EdHypoNeuro_mt10 = readRDS(url('https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/SeuratObj/EdHypoNeuro_mt10_integrated.rds'))

EdHypoNeuro_mt10@meta.data$region = 'tuberal'
EdHypoNeuro_mt10@meta.data$region[grep('HTHso',EdHypoNeuro_mt10@meta.data$sample)] = 'supraoptic'
EdHypoNeuro_mt10@meta.data$region[grep('HTHmn',EdHypoNeuro_mt10@meta.data$sample)] = 'mammillary'
EdHypoNeuro_mt10@meta.data$region[grep('HTHpo',EdHypoNeuro_mt10@meta.data$sample)] = 'preoptic'

## start hicat functions
ref.cl.df <- unique(data.frame(cluster_id = as.character(EdHypoNeuro_mt10@meta.data$integrated_snn_res.0.5) , cluster_label=paste('cluster',as.character(EdHypoNeuro_mt10@meta.data$integrated_snn_res.0.5),sep='_'), cluster_color = hue_pal()(38)[as.numeric(as.character(EdHypoNeuro_mt10@meta.data$integrated_snn_res.0.5) )+1], broad_type='Neuron'))

ref.cl.df <- arrange(ref.cl.df, cluster_id)
row.names(ref.cl.df) <- ref.cl.df$cluster_id

ref.cl <- setNames(factor(EdHypoNeuro_mt10@meta.data$region), rownames(EdHypoNeuro_mt10@meta.data))

EdHypoNeuro_mt10_cpm <- cpm(EdHypoNeuro_mt10@assays$RNA@counts)
norm.dat <- log2(EdHypoNeuro_mt10_cpm + 1)

norm.dat <- Matrix(norm.dat, sparse = TRUE)

de.param <- de_param(padj.th     = 0.01, 
                     lfc.th      = 1, 
                     low.th      = 1, 
                     q1.th       = 0.3,
                     q2.th       = NULL,
                     q.diff.th   = 0.3, 
                     de.score.th = 150,
                     min.cells = 10)

## skipping rm.eigen for now - didn't previously see batch effects

strict.param <- de_param(de.score.th = 500)

onestep.result <- onestep_clust(norm.dat, dim.method = "pca",  de.param = strict.param)

save(onestep.result,file='./Analysis/EdHypoNeuro_mt10_FirstIter_HiCat.rda',compress=TRUE)

iter.result <- iter_clust(norm.dat, dim.method = "pca",de.param = de.param, result = onestep.result)

save(iter.result,file='./Analysis/EdHypoNeuro_mt10_SecondIter.rda',compress=TRUE)

rd.dat <- t(norm.dat[iter.result$markers,])

merge.param <- de_param(de.score.th = 250) # The original value was 150. - rasing it should require more de genes between clusters.

merge.result <- merge_cl(norm.dat, 
                         cl = iter.result$cl, 
                         rd.dat = rd.dat,
                         de.param = merge.param) ## now 176 clusters 


save(merge.result,file='./Analysis/EdHypoNeuro_mt10_Merge_176clusters.rda',compress=TRUE)

## merge.result object saved here: 'https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/Analysis/EdHypoNeuro_mt10_Merge_176clusters.rda'
