library(Seurat)

flexsplit <- function(dat,char){
test=strsplit(as.character(dat),char,fixed=TRUE)
n.obs <- sapply(test, length)
seq.max <- seq_len(max(n.obs))
mat <- t(sapply(test, "[", i = seq.max))
mat
}

subplots <- function(sobj,colName,filename,raster=FALSE){
classes = unique(sobj@meta.data[colName])[,1]
pdf(filename)
for(i in classes){
sobj@meta.data[i]="Other"
sobj@meta.data[which(sobj@meta.data[colName]==i),i]=i
print(DimPlot(sobj, reduction = "umap", group.by = i,order=c(i,'other'),raster=raster) + ggplot2::labs(title=i))
}
dev.off()
}

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

MM2HSref = read.csv('https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/Reference_Files/MM2HS_EG100.csv')
colnames(MM2HSref) = c('MM_ID','MM_Gene','HS_ID','HS_Gene')
HS2MMref = read.csv('https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/Reference_Files/HS2MM_EG100.csv')
colnames(HS2MMref) = c('HS_ID','HS_Gene','MM_ID','MM_Gene')


HStoMM <- function(x){
    tmpind = match(as.character(x),HS2MMref$HS_Gene)
return(as.character(HS2MMref$MM_Gene[tmpind]))
}

MMtoHS <- function(x){
    tmpind = match(as.character(x),MM2HSref$MM_Gene)
return(as.character(MM2HSref$HS_Gene[tmpind]))
}

HSTF = read.csv('https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/Reference_Files/Human_TF.csv')

NPlist = read.csv('https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/Reference_Files/Neuropeptide_list.csv',header=FALSE) ## recieved from Hannah on 7/15/20

NPlist = MMtoHS(as.character(NPlist[,1]))


TF_NP = read.csv('https://data.nemoarchive.org/biccn/grant/u01_devhu/kriegstein/transcriptome/scell/10x_v2/human/processed/analysis/Herb_2022_Hypothalamus/Reference_Files/TF_NP_list_from_Moffitt.csv') #direct from Moffitt suppmental materials 

for(i in 1:ncol(TF_NP)){
    TF_NP[,i] = MMtoHS(TF_NP[,i])
} ## human convention 


NPlist2 = unique(c(na.omit(as.character(TF_NP$Neuropeptides)),NPlist))

getDescendants<-function(tree,node,curr=NULL){
  if(is.null(curr)) curr<-vector()
  daughters<-tree$edge[which(tree$edge[,1]==node),2]
  curr<-c(curr,daughters)
  w<-which(daughters>=length(tree$tip))
  if(length(w)>0) for(i in 1:length(w))
    curr<-getDescendants(tree,daughters[w[i]],curr)
  return(curr)
}

getTips = function(tree=NA,node=NA){
tmp = getDescendants(tree,node)
tmp2 = tmp[which(tmp<=length(tree$tip.label))]
return(tmp2)
}
