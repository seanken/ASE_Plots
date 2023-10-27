library(scAlleleExpression)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(purrr)
library(stringi)
library(dplyr)
library(tidyr)
library(DropletUtils)

qc=readRDS("../qc.picked.genes.RDS")
genes=qc[qc[,"Picked"],"ID"]

lst=system("ls /stanley/levin_asap_storage/6*/n*21*/CellRanger/*hybsel/s1/outs/m*h5",intern=T)
mol.info=read10xMolInfo(lst[1])
dat=mol.info[[1]]
gene_names=mol.info[[2]]

head(dat)
dat["In"]=gene_names[dat[,"gene"]] %in% genes
dat=dat[sample(1:dim(dat)[1],100000),]

dat1=dat

mol.info=read10xMolInfo(lst[2])
dat=mol.info[[1]]
gene_names=mol.info[[2]]

head(dat)
dat["In"]=gene_names[dat[,"gene"]] %in% genes
dat=dat[sample(1:dim(dat)[1],100000),]
dat2=dat

dat=rbind(dat1,dat2)
dat=data.frame(dat)

dat["Type"]="In Chosen Gene"
dat[dat[,"In"]==F,"Type"]="Not In Chosen Gene"

p=ggplot(dat,aes(x=Type,y=reads,fill=In))+geom_violin(scale="width")+scale_y_log10()+xlab("")+ylab("Reads per UMI")+theme(legend.position="none")
ggsave("Reads.per.UMI.pdf",p)