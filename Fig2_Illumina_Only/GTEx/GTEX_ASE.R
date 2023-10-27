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

##Add in code for sig AI genes
ASE=apply(fils,1,function(x){
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

names(ASE)=fils[,"samp"]
saveRDS(ASE,"ASE.RDS")
bulk=readRDS('GTEx.bulk.RDS')
for(i in names(ASE)){ASE[[i]]["Samp"]=i}
dat=do.call(rbind,ASE)
dat=dat[,c("Gene","Samp","padj","Estimate")]
colnames(dat)[1]="SYMBOL"
bulk<-bulk %>% gather(Col,Val,YJ89_ref,YJ89_alt,Z93S_ref,Z93S_alt,YJ89_binom_padj,Z93S_binom_padj)
bulk["Samp"]=map_chr(bulk[,"Col"],function(x) strsplit(x,"_")[[1]][1])
bulk["valtype"]=map_chr(bulk[,"Col"],function(x) strsplit(x,"_")[[1]][2])
bulk=bulk[,c("SYMBOL","Val","Samp","valtype")]
bulk <-bulk %>% spread(valtype,Val,fill=0)
bulk=bulk[bulk$binom!=-1,]
comb=inner_join(dat,bulk)
comb["Bulk"]=log((comb[,"ref"]+1)/(comb[,"alt"]+1),2)
comb=comb[!is.na(comb$padj),]
#p=ggplot(comb[comb$padj<.05,],aes(x=Estimate,y=Bulk,color=Samp))+geom_point()+geom_hline(yintercept=0,linetype="dotted")+geom_vline(xintercept=0,linetype="dotted")+geom_abline(slope=1,linetype="dotted")

p=ggplot(comb[comb$padj<.05,],aes(x=Estimate,y=Bulk,color=Samp))+geom_point()+geom_hline(yintercept=0,linetype="dotted")+geom_vline(xintercept=0,linetype="dotted")+geom_abline(slope=1,linetype="dotted")+ylab("Effect Size From Bulk")+xlab("Effect Size From Single Nuclei")
ggsave("CompareToBulk.pdf",p)
saveRDS(comb,"comb.RDS")

