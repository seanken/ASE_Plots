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

print("Get reads")
numReads=map_dbl(dirs,function(x){
    tab=read.csv(paste(x,"/STARSolo/output/resultsSolo.out/GeneFull/Summary.csv",sep=""),header=F)
    tab[1,2]
})

print("Get files")
fils=map(dirs,function(x){
    fil=map_chr(1:9,function(y) paste(x,"/Downsample/counts.",y,".txt",sep=""))
    fil=c(fil,paste(x,"/AlleleCounts/counts.txt",sep=""))
    DS=.1*1:10
    dat=data.frame(fil,DS)
})

for(i in names(fils)){
    fils[[i]]["samp"]=strsplit(i,"_")[[1]][1]
    fils[[i]]["hybType"]=strsplit(i,"_")[[1]][2]
}


fils=do.call(rbind,fils)

fils["num"]=map_chr(fils[,"fil"],function(x) system(paste("head",x,"| wc -l"),intern=T))
fils2=fils
fils=fils[fils$num=="10",] ##Will need to fix

fils["nPhased"]=map_int(fils[,"fil"],function(x){
    
    tab=read.table(x)
    print(head(tab))
    tab=tab[tab[,3]!="Ambig",]
    tab=tab[tab[,2] %in% genes,]
    tot=sum(tab[,4])
    print(tot)
    return(tot)
})

fils<-fils %>% unite(Nam,samp,hybType,sep="_",remove=F)
fils["NumReads"]=fils[,"DS"]*numReads[fils[,"Nam"]]
fils[fils$hybType=="hyb","hybType"]="Selected"
fils[fils$hybType=="nohyb","hybType"]="Not Selected"
p=ggplot(fils,aes(x=NumReads,y=nPhased,color=hybType))+geom_point()+geom_line()+ylab("Number phased UMI in Target Genes")+xlab("Number Sequenced Reads (Illumina)")+theme(legend.title=element_blank())+facet_wrap(~samp,ncol=2)
ggsave("Downsample.Selection.Reads.PhasedUMI.pdf",p,width=14)

##Add in code for sig AI genes
fils["nSig"]=apply(fils,1,function(x){
    print(x["DS"])
    tab=read.table(x["fil"])
    tab["Name"]=sub("^",paste(x["samp"],"_",sep=""),tab[,1])
    print(dim(tab))
    tab=tab[tab$Name %in% names(seur@active.ident),]
    print(dim(tab))
    colnames(tab)[2]="Gene"
    tab["Sample"]="samp"
    colnames(tab)[3]="Allele"
    colnames(tab)[4]="Count"
    tab["SNP"]=tab["Gene"]
    tab=tab[tab[,2] %in% genes,]
    tab[tab$Allele=="All1","Allele"]="alt"
    tab[tab$Allele=="All2","Allele"]="ref"
    num<-tab %>% group_by(Gene) %>% summarise(numAlt=sum(Allele=="alt"),numRef=sum(Allele=="ref")) %>% as.data.frame()
    tab=tab[tab$Gene %in% num[num$numAlt>0 & num$numRef>0,"Gene"],]
    out=TestSNP_aod(tab,minCount=10,minSamp=0)
    return(sum(out$padj<.05,na.rm=T))
})



p=ggplot(fils,aes(x=NumReads,y=nSig,color=samp))+geom_point()+geom_line()+ylab("Number Genes with Significant AI")+xlab("Number Sequenced Reads (Illumina)")+theme(legend.title=element_blank())
ggsave("Downsample.Reads.SigGenes.pdf",p)