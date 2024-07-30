library(scAlleleExpression)
library(Seurat)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(purrr)
library(tidyr)
library(dplyr)

print("Load Seurat Data")
seur=readRDS("../../Fig2_Illumina_Only/UMAP/seur.RDS")

qc=readRDS("../qc.picked.genes.RDS")
genes=qc[qc[,"Picked"],"Gene"]

dirs=system("ls /stanley/levin_asap_storage/612-eqtl/DS_ForPaper/GTExData/*Hyb*/output | grep : | sed 's/://g'",intern=T)
names(dirs)=c("YJ89_hyb","YJ89_nohyb","Z93S_hyb","Z93S_nohyb")
qc=readRDS("../qc.picked.genes.RDS")
genes=qc[qc[,"Picked"],"Gene"]

out=do.call(rbind,map(names(dirs),function(y){
    print(y)
    x=paste(dirs[y],"/AlleleCounts/counts.txt",sep="")
    tab=read.table(x)
    print(head(tab))
    tab=tab[tab[,3]!="Ambig",]
    tab=tab[tab[,2] %in% genes,]
    tab<-tab %>% group_by(V2) %>% summarise(Tot=sum(V4)) %>% as.data.frame()
    tab["Test"]=y
    return(tab)
}))

out["Sample"]=map_chr(out[,"Test"],function(x) strsplit(x,"_")[[1]][1])
out["Type"]=map_chr(out[,"Test"],function(x) strsplit(x,"_")[[1]][2])

out=out[,c("V2","Tot","Sample","Type")]
out<-out %>% spread(Type,Tot,fill=0) %>% as.data.frame()
saveRDS(out,"out.for.bygene.RDS")
p= ggplot(out,aes(x=hyb,y=nohyb,color=Sample))+geom_point()+geom_abline(slope=1,intercept=0,linetype="dotted")+coord_flip()+xlab("Number Phased UMI With Selection")+ylab("Number Phased UMI No Selection")+scale_x_log10()+scale_y_log10()
ggsave("PerGene.pd",p)
