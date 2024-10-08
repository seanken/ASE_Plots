source("/stanley/levin_dr/ssimmons/SingleCell3/load_Seurat.R")
library(ggplot2)
library(qs)
#dat=readRDS("/stanley/levin_asap_storage/612-eqtl/SingleNuc_data/UpdatedAIData/RunNewData/DataFromClemens/TransEQTL/expression.1perc.RDS")
print("Load data")
dat=qread("/stanley/levin_asap/ssimmons/eQTL/PlotsForPaper/Fig5_SNP_Level_PD/UMAP/dat.expression.qs")
meta=qread("/stanley/levin_asap_storage/ssimmons/ClemensUpdated/meta.with.ASEInfo.qs")
#meta=readRDS("/stanley/levin_asap_storage/612-eqtl/SingleNuc_data/UpdatedAIData/RunNewData/DataFromClemens/TransEQTL/meta.RDS")
print("Make seurat")
seur=dir10X(dat=dat,minGenes=0,regress=c())
print("Run UMAP")
seur=RunUMAP(seur,dims=1:20)
seur@meta.data["CellType"]=meta[,"MajorCellTypes"]
print("Make plot")
p=UMAPPlot(seur,label=T,group.by="CellType")
print("SAVE")
qsave(p,"UMAP.qs")
ggsave("UMAP.pdf",p,width=14)
qsave(seur,"seur.qs")
