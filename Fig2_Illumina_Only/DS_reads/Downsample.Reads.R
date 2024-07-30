library(scAlleleExpression)
library(Seurat)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(purrr)

print("Load Seurat Data")
seur=readRDS("../UMAP/seur.RDS")
dirs=c("/stanley/levin_asap_storage/612-eqtl/DS_ForPaper/GTExData/samp_YJ89_130/output","/stanley/levin_asap_storage/612-eqtl/DS_ForPaper/GTExData/samp_Z93S_130/output")
names(dirs)=c("YJ89","Z93S")


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
    fils[[i]]["samp"]=i
}

fils=do.call(rbind,fils)

fils["nPhased"]=map_int(fils[,"fil"],function(x){
    tab=read.table(x)
    print(head(tab))
    tab=tab[tab[,3]!="Ambig",]
    return(sum(tab[,4]))
})

fils["NumReads"]=fils[,"DS"]*numReads[fils[,"samp"]]
p=ggplot(fils,aes(x=NumReads,y=nPhased,color=samp))+geom_point()+geom_line()+ylab("Number phased UMI")+xlab("Number Sequenced Reads (Illumina)")+theme(legend.title=element_blank())
ggsave("Downsample.Reads.PhasedUMI.pdf",p)

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
    tab[tab$Allele=="All1","Allele"]="alt"
    tab[tab$Allele=="All2","Allele"]="ref"
    num<-tab %>% group_by(Gene) %>% summarise(numAlt=sum(Allele=="alt"),numRef=sum(Allele=="ref")) %>% as.data.frame()
    tab=tab[tab$Gene %in% num[num$numAlt>0 & num$numRef>0,"Gene"],]
    out=TestSNP_aod(tab,minCount=10,minSamp=0)
    return(sum(out$padj<.05,na.rm=T))
})


saveRDS(fils,"fils.RDS")
p=ggplot(fils,aes(x=NumReads,y=nSig,color=samp))+geom_point()+geom_line()+ylab("Number Genes with Significant AI")+xlab("Number Sequenced Reads (Illumina)")+theme(legend.title=element_blank())
ggsave("Downsample.Reads.SigGenes.pdf",p)
