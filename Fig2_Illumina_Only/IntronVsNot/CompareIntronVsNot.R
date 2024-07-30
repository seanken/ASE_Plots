library(scAlleleExpression)
library(Seurat)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(purrr)

print("Load Seurat Data")
seur=readRDS("../../Fig2_Illumina_Only/UMAP/seur.RDS")

print("Set up data frame")
lst=system("ls /stanley/levin_asap_storage/612-eqtl/DS_ForPaper/GTExData/samp_* | grep -v Hyb | grep : | sed 's/://g' | grep '0$\\|5$' | grep 130",intern=T)
counts=sub("$","/output/AlleleCounts/counts.txt",lst)
len=map_dbl(lst,function(x){s=strsplit(x,"_")[[1]];as.numeric(s[length(s)])})
samp=map_chr(lst,function(x){s=strsplit(x,"_")[[1]];s[length(s)-1]})

fils=data.frame("fil"=counts,"len"=len,"samp"=samp)
fils["QC"]=counts=sub("$","/output/STARSolo/output/resultsSolo.out/GeneFull/Summary.csv",lst)



fils["TotUMI"]=map_dbl(fils[,"QC"],function(x){read.csv(x,header=F)[15,2]})
#fils["percPhase"]=100*fils[,"nPhased"]/fils[,"TotUMI"]

print("Perform ASE testing Illumina")
ASE_ill=apply(fils,1,function(x){
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
    return(out)
})

#seur=readRDS("../DSReads/seur.MAS.RDS")

fils["fil_MAS"]=c("/stanley/levin_asap_storage/612-eqtl/DS_ForPaper/GTExData/NoIntrons/samp_nointron_YJ89_130/output/AlleleCounts/counts.txt","/stanley/levin_asap_storage/612-eqtl/DS_ForPaper/GTExData/NoIntrons/samp_nointron_Z93S_130/output/AlleleCounts/counts.txt")

print("Perform ASE testing MAS-Seq")
ASE_MAS=apply(fils,1,function(x){
    print(x["len"])
    tab=read.table(x["fil_MAS"])
    tab["Name"]=sub("^",paste(x["samp"],"_",sep=""),tab[,1])
    print(dim(tab))
    tab=tab[tab$Name %in% names(seur@active.ident),]
    print(dim(tab))
    print(length(unique(tab$Name)))
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
    return(out)
})

comb1=inner_join(ASE_ill[[1]],ASE_MAS[[1]],by="Gene")
comb2=inner_join(ASE_ill[[2]],ASE_MAS[[2]],by="Gene")

comb1["Name"]=fils[1,"samp"]
comb2["Name"]=fils[2,"samp"]

comb=rbind(comb1,comb2)
saveRDS(comb,"comb.RDS")
comb=comb[!is.na(comb$padj.x) & !is.na(comb$padj.y),]
p=ggplot(comb[comb$padj.x<.05 | comb$padj.y<.05,],aes(x=Estimate.x,y=Estimate.y,color=Name))+geom_point()+xlab("Effect size with introns")+ylab("Effect size no introns")+geom_vline(xintercept=0,linetype="dotted")+geom_hline(yintercept=0,linetype="dotted")
ggsave("Compare.pdf",p)
