
library(scAlleleExpression)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

ref=NULL
if(!file.exists("ref.RDS"))
{
    print("Prep ref data")
    meta=readRDS("/stanley/levin_asap_storage/612-eqtl/SingleNuc_data/UpdatedAIData/RunNewData/DataFromClemens/meta.with.files.RDS")
    set.seed(1)
    meta=meta[sample(rownames(meta),50000),]
    meta["Raw"]=sub("$","/output/STARSolo/output/resultsSolo.out/GeneFull/raw",meta[,"dir_ASE"])
    dat=LoadExpression(meta, cbc_col = "CBC", raw_col = "Raw",samp_col="batchSamp")
    source("/stanley/levin_dr/ssimmons/SingleCell3/load_Seurat.R")
    seur_ref=dir10X(dat=dat,minGenes=0)
    seur_ref=RunUMAP(seur_ref,dims=1:20,return.model=T)
    seur_ref@meta.data["CellType"]=meta[,"MajorCellTypes"]
    source("/stanley/levin_diamond/UpdatedTerra/Results/Azimuth.R")
    ref=GenerateReference(seur_ref,"CellType")
    rm(seur_ref)
    saveRDS(dat,"dat.for.Clemens.RDS")
    saveRDS(ref,"ref.RDS")
}
else{
    ref=readRDS("ref.RDS")
}


print("Load data")
fils=c("/stanley/levin_asap_storage/612-eqtl/DS_ForPaper/GTExData/samp_YJ89_130/output/STARSolo/output/resultsSolo.out/GeneFull/filtered","/stanley/levin_asap_storage/612-eqtl/DS_ForPaper/GTExData/samp_Z93S_130/output/STARSolo/output/resultsSolo.out/GeneFull/filtered")
names(fils)=c("YJ89","Z93S")
dat=ReadInSTARSolo(fils)

source("/stanley/levin_dr/ssimmons/SingleCell3/SharedVariable.v2.R")
print("Cluster")
seur<-CreateSeuratObject(dat,"Seurat",min.features=250)
seur<-NormalizeData(seur,normalization.method="LogNormalize",scale.factor=1000000)
seur<-SharedVariable(seur,x.low.cutoff=1,x.high.cutoff=5,minNum=0,batch="orig.ident",minCells=10)
#seur<-FindVariableFeatures(seur)
seur<-ScaleData(seur,features=seur@assays$RNA@var.features,vars.to.regress=c("nFeature_RNA","orig.ident"))
seur<-RunPCA(seur,npcs=60)
#library(harmony)
#seur=RunHarmony(seur,"orig.ident",dims=1:20)
library(SeuratWrappers)
seur<- RunFastMNN(object.list = SplitObject(seur, split.by = "orig.ident"))
seur=FindNeighbors(seur,dims=1:20,reduction="mnn")
seur=FindClusters(seur)
seur=RunUMAP(seur,dims=1:20,reduction="mnn")
print("Next")
fils=sub("/STARSolo/output/resultsSolo.out/GeneFull/filtered","",fils,fixed=T)
meta=LoadPipeline(seur@meta.data,fils)
seur@meta.data=meta

print("Get doublets")
source("/stanley/levin_dr/ssimmons/SingleCell3/RemoveDoublet.R")
seur=DoubletScores(seur)

print("Markers")
genes=c("GAD1","SLC17A7","SNAP25","PLP1","PDGFRA","SLC1A3","CSF1R","scds","FLT1")
p=FeaturePlot(seur,genes,ncol=3)
ggsave("Markers.pdf",p,width=14,height=10)


print("Label Cell Types")
seur@meta.data["Cluster"]=as.numeric(as.character(seur@active.ident))+1
lst=rep("Excitatory",30)
lst[c(29,22,19,28,20,21,25,30
)]="Doublet"
lst[c(11,12,14,18)]="Inhibitory"
lst[c(3,9)]="ODC"
lst[6]="OPC"
lst[5]="Microglia"
lst[c(2,15)]="Astrocyte"
lst[23]="Vascular" ##EPAS1+
seur@meta.data["CellType"]=lst[as.numeric(seur@meta.data$Cluster)]


#print("Run Azimuth")
#source("/stanley/levin_diamond/UpdatedTerra/Results/Azimuth.R")
#meta=RunAzimuth(seur,ref)



print("Simple plots")
p=UMAPPlot(seur,label=T,group.by="orig.ident")
ggsave("UMAP.sample.pdf",p,width=14,height=10)

p=UMAPPlot(seur,label=T,group.by="CellType")
ggsave("UMAP.celltype.pdf",p,width=14,height=10)

print("Save")
saveRDS(seur,"seur.RDS")
seur=subset(seur,CellType!="Doublet")
saveRDS(seur,"seur.no.doublets.RDS")
print("Done")




