
library(scAlleleExpression)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(qs)
source("/stanley/levin_dr/ssimmons/SingleCell3/Run_scGBM.R")
library(Seurat)



dat=qread("counts.with.cells.with.Revio.qs")
##used out2=map(out,function(dat){tab=emptyDropsCellRanger(dat,n.expected.cells=10000,umi.min=250);cells=rownames(tab)[tab$FDR <= .01 & !is.na(tab[,"FDR"]) ];dat=dat[,cells]})
names(dat)=rev(c("YJ89","Z93S"))
for(i in names(dat)){colnames(dat[[i]])=sub("^",paste(i,"_",sep=""),colnames(dat[[i]]))}
inter=intersect(rownames(dat[[1]]),rownames(dat[[2]]))
mat=cbind(dat[[1]][inter,],dat[[2]][inter,])


source("/stanley/levin_dr/ssimmons/SingleCell3/SharedVariable.v2.R")
print("Cluster")
seur<-CreateSeuratObject(mat,"Seurat",min.features=150)
seur<-NormalizeData(seur,normalization.method="LogNormalize",scale.factor=10000)
seur<-SharedVariable(seur,x.low.cutoff=.1,x.high.cutoff=.5,minNum=0,batch="orig.ident",minCells=10)
#seur=Run_scGBM(seur,numDims=20,batch="orig.ident")
dat=seur[[assay]]@counts[seur[[assay]]@var.features,]
dat=as.matrix(dat)
out.proj <- gbm.sc(dat,M=numDims,subset=10000)
seur[["gbm"]] <- CreateDimReducObject(embeddings=out.proj$V,key="GBM_")
rownames(seur@reductions$gbm@cell.embeddings)=names(seur@active.ident)
seur=RunUMAP(seur,dims=1:numDims,reduction="gbm")
seur=FindNeighbors(seur,dims=1:numDims,reduction="gbm");seur=FindClusters(seur)

source("/stanley/levin_dr/ssimmons/SingleCell3/RemoveDoublet.R")
seur=DoubletScores(seur)
qsave(seur,"seur.MAS.qs")
seur=subset(seur,seurat_clusters!=11 & seurat_clusters!=12) #remove doiublets
seur=subset(seur,scds<1)
qsave(seur,"seur.nodoub.MAS.qs")

print("Markers")
genes=c("GAD1","SLC17A7","RBFOX3","PLP1","PDGFRA","SLC1A3","CSF1R","scds","FLT1")
p=FeaturePlot(seur,genes,ncol=3)
ggsave("Markers.scGBM.MAS.pdf",p,width=14,height=10)

p=UMAPPlot(seur,label=T,group.by="orig.ident")
ggsave("UMAP.sample.scGBM.MAS.pdf",p,width=14,height=10)

seur@meta.data["CellType"]="Excitatory"
seur@meta.data[seur@meta.data[,"seurat_clusters"] %in% c(7,10),"CellType"]="Inhibitory"
seur@meta.data[seur@meta.data[,"seurat_clusters"]==8,"CellType"]="OPC"
seur@meta.data[seur@meta.data[,"seurat_clusters"]==0,"CellType"]="ODC"
seur@meta.data[seur@meta.data[,"seurat_clusters"]==6,"CellType"]="Microglia"
seur@meta.data[seur@meta.data[,"seurat_clusters"]==2,"CellType"]="Astrocytes"
seur@meta.data[seur@meta.data[,"seurat_clusters"]==13,"CellType"]="Vascular"

p=UMAPPlot(seur,label=T,group.by="CellType")
ggsave("CellType.MAS.pdf",p)
qsave(seur,"seur.nodoub.MAS.qs")