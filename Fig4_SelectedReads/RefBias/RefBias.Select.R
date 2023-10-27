library(scAlleleExpression)
library(Seurat)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(purrr)
library(dplyr)
library(tidyr)

dirs=system("ls /stanley/levin_asap_storage/612-eqtl/DS_ForPaper/GTExData/*Hyb*/output | grep : | sed 's/://g'",intern=T)
names(dirs)=c("YJ89_hyb","YJ89_nohyb","Z93S_hyb","Z93S_nohyb")
dirs=dirs[c("YJ89_hyb","Z93S_hyb")]
qc=readRDS("../qc.picked.genes.RDS")
genes=qc[qc[,"Picked"],"Gene"]

out=map(dirs,function(x){
    print(x)
    dat=read.table(paste(x,"/SNPLevelAlleleCounts/counts.txt",sep=""))
    dat<-dat %>% group_by(V2,V3) %>% summarise(Count=sum(V4)) %>% spread(V3,Count,fill=0) %>% as.data.frame()
    colnames(dat)[1]="SNP_short"
    tab=read.table(paste(x,"/SNPLevelCounts/snps.bed",sep=""))
    tab2=read.table(paste(x,"/SNPLevelCounts/comb.bed",sep=""))
    tab2=tab2[tab2[,2] %in% genes,]
    snps=unique(tab2[,1])
    tab=tab[tab[,"V4"] %in% snps,]
    tab<-tab %>% unite(SNP_short,V1,V2,sep=":")
    comb=inner_join(dat,tab,by="SNP_short")
    comb=comb[comb$V5!=".|.",]
    comb["ref"]=comb[,"All1"]
    comb["alt"]=comb[,"All2"]
    comb[comb$V5=="1|0","ref"]=comb[comb$V5=="1|0","All1"]
    comb[comb$V5=="1|0","alt"]=comb[comb$V5=="1|0","All2"]
    comb["rat"]=comb["ref"]/(comb[,"ref"]+comb[,"alt"])
    comb["tot"]=comb[,"ref"]+comb[,"alt"]
    comb=comb[comb$tot>10,]
    return(comb)
})

for(i in names(out)){
    out[[i]]["Name"]=i
    out[[i]]["samp"]=strsplit(i,"_")[[1]][1]
    out[[i]]["hybType"]=strsplit(i,"_")[[1]][1]
}

dat=do.call(rbind,out)
saveRDS(dat,"ratios.RDS")

p=ggplot(dat,aes(x=samp,y=100*rat,fill=samp))+geom_violin(scale="width")+ylab("Percent Phased UMI from Reference")+xlab("Sample")+theme(legend.position="none")
ggsave("RefBias.Select.pdf",p)
#dat=dat[dat$rat==0 | dat$rat==1,]
#tab<-dat %>% group_by(samp,hybType,rat) %>% summarise(Num=length(samp)) %>% as.data.frame()
#tab["Type"]="Reference only"
##tab[tab$rat==0,"Type"]="Alternative only"
#p=ggplot(tab,aes(x=Type,fill=samp,y=Num))+geom_bar(stat="identity",position="dodge")+xlab("")+ylab("Number of SNPs")+facet_wrap(~)
#ggsave("Number.Ref.vs.Alt.pdf")