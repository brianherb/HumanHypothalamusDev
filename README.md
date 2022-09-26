# Single-cell genomics reveals region-specific developmental trajectories underlying neuronal diversity in the prenatal human hypothalamus

## Code in this repository emcompases the analysis performed for the manuscript "Single-cell genomics reveals region-specific developmental trajectories underlying neuronal diversity in the prenatal human hypothalamus"

The scripts follow the order presented in the manuscript:

environment_droplet.yml - Conda yml file - can be used to recreate conda environment used in this analysis

HypoPub_functions_reference.R - Collection of functions and reference files used in the following scripts

Read_in_Data.R - converts MTX counts files into Seurat objects

Integrate_Samples.R - Individual Seurat objects were integrated together for various downstream analyses

HiCat_Clustering_HumanAdult.R - HiCat clustering of Adult Human Neurons. 

Fig1_HumanEmbryonic.R - Identification of cell types in Human Embryoic samples 

Fig2_HumanNeuron.R - Construction of developmental lineages across embryonic and adult hypothalamic neurons 

Fig3_HumanMouseNeuron.R - Comparison of mouse and human hypothalamic neurons. 

Fig4_HumanHypoCtxGE.R - Comparison of transcription factor expression and gene regulatory networks across hypothalamus, ganglionic eminence and cortex in the developing human brain
