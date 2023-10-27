library(scAlleleExpression)
library(Seurat)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(purrr)

print("Load Seurat Data")
seur=readRDS("../UMAP/seur.RDS")

print("Set up data frame")
lst=system("ls /stanley/levin_asap_storage/612-eqtl/DS_ForPaper/GTExData/samp_* | grep -v Hyb | grep : | sed 's/://g' | grep '0$\\|5$' ",intern=T)
counts=sub("$","/output/AlleleCounts/counts.txt",lst)
len=map_dbl(lst,function(x){s=strsplit(x,"_")[[1]];as.numeric(s[length(s)])})
samp=map_chr(lst,function(x){s=strsplit(x,"_")[[1]];s[length(s)-1]})

fils=data.frame("fil"=counts,"len"=len,"samp"=samp)
fils["QC"]=counts=sub("$","/output/STARSolo/output/resultsSolo.out/GeneFull/Summary.csv",lst)

fils["nPhased"]=map_int(fils[,"fil"],function(x){
    tab=read.table(x)
    print(head(tab))
    tab=tab[tab[,3]!="Ambig",]
    return(sum(tab[,4]))
})

fils["TotUMI"]=map_dbl(fils[,"QC"],function(x){read.csv(x,header=F)[15,2]})
fils["percPhase"]=100*fils[,"nPhased"]/fils[,"TotUMI"]

##Add in code for sig AI genes
fils["nSig"]=apply(fils,1,function(x){
    print(x["len"])
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
fils["len"]=factor(fils[,"len"],c(55,80,130,180))
p=ggplot(fils,aes(x=len,y=percPhase,fill=samp))+geom_bar(position="dodge",stat="identity")+ylab("Percent of UMI that are Phased")+xlab("Length of Read 2 (bp)")+theme(legend.title=element_blank())
ggsave("Trim.Reads.PercentPhasedUMI.pdf",p)


p=ggplot(fils,aes(x=len,y=nSig,fill=samp))+geom_bar(position="dodge",stat="identity")+ylab("Number Genes with Sig AI")+xlab("Length of Read 2 (bp)")+theme(legend.title=element_blank())
ggsave("Trim.Reads.SigGenes.pdf",p)