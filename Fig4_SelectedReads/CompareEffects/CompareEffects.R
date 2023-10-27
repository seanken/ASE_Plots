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

dirs=system("ls /stanley/levin_asap_storage/612-eqtl/DS_ForPaper/GTExData/*Hyb*/output | grep : | sed 's/://g'",intern=T)
names(dirs)=c("YJ89_hyb","YJ89_nohyb","Z93S_hyb","Z93S_nohyb")
#dirs=dirs[c("YJ89_hyb","Z93S_hyb")]
qc=readRDS("../qc.picked.genes.RDS")
genes=qc[qc[,"Picked"],"Gene"]

out=map(names(dirs),function(x){
    print(x)
    tab=read.table(paste(dirs[x],"/AlleleCounts/counts.txt",sep=""))
    print(head(tab))
    tab=tab[tab[,3]!="Ambig",]
    tab=tab[tab[,2] %in% genes,]
    samp=strsplit(x,"_")[[1]][1]
    tab["Name"]=sub("^",paste(samp,"_",sep=""),tab[,1])
    print(dim(tab))
    tab=tab[tab$Name %in% names(seur@active.ident),]
    print(dim(tab))
    colnames(tab)[2]="Gene"
    tab["Sample"]="samp"
    colnames(tab)[3]="Allele"
    colnames(tab)[4]="Count"
    tab["SNP"]=tab["Gene"]
    #tab=tab[tab[,2] %in% genes,]
    tab[tab$Allele=="All1","Allele"]="alt"
    tab[tab$Allele=="All2","Allele"]="ref"
    num<-tab %>% group_by(Gene) %>% summarise(numAlt=sum(Allele=="alt"),numRef=sum(Allele=="ref")) %>% as.data.frame()
    tab=tab[tab$Gene %in% num[num$numAlt>0 & num$numRef>0,"Gene"],]
    mrk=TestSNP_aod(tab,minCount=10,minSamp=0)
    print(head(mrk))
    return(mrk)
})

print("Prep")
names(out)=names(dirs)
comb1=inner_join(out[["YJ89_hyb"]],out[["YJ89_nohyb"]],by="Gene")
comb2=inner_join(out[["Z93S_hyb"]],out[["Z93S_nohyb"]],by="Gene")
comb1["Sample"]="YJ89"
comb2["Sample"]="Z93S"
comb=rbind(comb1,comb2)
saveRDS(comb,"AI.RDS")


comb["Sig"]="Not"
comb[comb$padj.x<.05 & !is.na(comb$padj.x),"Sig"]="FDR<.05, Selected"
comb[comb$padj.y<.05 & !is.na(comb$padj.y),"Sig"]="FDR<.05, Not Selected"
comb[comb$padj.x<.05 & comb$padj.y<.05,"Sig"]="FDR<.05, Both"
print("plot")
p=ggplot(comb,aes(x=Estimate.x,y=Estimate.y,color=Sig))+geom_point()+facet_wrap(~Sample)+ylab("Estimated AI Effect Size, No Selection")+xlab("Estimated AI Effect Size, No Selection")+geom_abline(slope=1,linetype="dotted")
ggsave("CompareEffect.pdf",p,width=14)

dat=comb
dat["has_eQTL"]="No"
dat[dat$Gene %in% qc[qc$Reason_Picked=="QTL" & qc$Picked,"Gene"],"has_eQTL"]="Yes"
dat=dat %>% group_by(has_eQTL,Sample) %>% summarise(PercentSig=100*mean(padj.x<.05 & !is.na(padj.x)),Num=length(padj.x)) %>% as.data.frame()
p=ggplot(dat,aes(x=has_eQTL,y=PercentSig,fill=Sample))+geom_bar(stat="identity",position="dodge")+ylab("Percent Genes with Significant AI")+xlab("Gene has known eQTL in it")
ggsave("BarPlot.eQTL.pdf")
