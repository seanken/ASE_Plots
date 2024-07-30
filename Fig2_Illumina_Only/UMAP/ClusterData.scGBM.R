
library(scAlleleExpression)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(qs)
source("/stanley/levin_dr/ssimmons/SingleCell3/Run_scGBM.R")


print("Load data")
fils=c("/stanley/levin_asap_storage/612-eqtl/DS_ForPaper/GTExData/samp_YJ89_130/output/STARSolo/output/resultsSolo.out/GeneFull/filtered","/stanley/levin_asap_storage/612-eqtl/DS_ForPaper/GTExData/samp_Z93S_130/output/STARSolo/output/resultsSolo.out/GeneFull/filtered")
names(fils)=c("YJ89","Z93S")
dat=ReadInSTARSolo(fils)

source("/stanley/levin_dr/ssimmons/SingleCell3/SharedVariable.v2.R")
print("Cluster")
seur<-CreateSeuratObject(dat,"Seurat",min.features=250)
seur<-NormalizeData(seur,normalization.method="LogNormalize",scale.factor=10000)
seur<-SharedVariable(seur,x.low.cutoff=.1,x.high.cutoff=.5,minNum=0,batch="orig.ident",minCells=10)
seur=Run_scGBM(seur,numDims=20,batch="orig.ident")


print("Next")
fils=sub("/STARSolo/output/resultsSolo.out/GeneFull/filtered","",fils,fixed=T)
meta=LoadPipeline(seur@meta.data,fils)
seur@meta.data=meta

print("Get doublets")
source("/stanley/levin_dr/ssimmons/SingleCell3/RemoveDoublet.R")
seur=DoubletScores(seur)

qsave(seur,"seur.GBM.qs")

print("Markers")
genes=c("GAD1","SLC17A7","SNAP25","PLP1","PDGFRA","SLC1A3","CSF1R","scds","FLT1")
p=FeaturePlot(seur,genes,ncol=3)
ggsave("Markers.scGBM.pdf",p,width=14,height=10)

p=UMAPPlot(seur,label=T,group.by="orig.ident")
ggsave("UMAP.sample.scGBM.pdf",p,width=14,height=10)


tab=qread("CT.qs")

seur@meta.data["CellType"]=tab[as.numeric(as.character(seur@meta.data[,"seurat_clusters"]))+1,"CellType"]

p=UMAPPlot(seur,label=T,group.by="CellType")
ggsave("UMAP.sample.celltype.pdf",p,width=14,height=10)
qsave(seur,"seur.GBM.qs")
seur=subset(seur,CellType!="Doublet")
qsave(seur,"seur.GBM.qs")
#print("Run Azimuth")
#source("/stanley/levin_diamond/UpdatedTerra/Results/Azimuth.R")
#meta=RunAzimuth(seur,ref)

#



