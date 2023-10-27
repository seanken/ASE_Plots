library(scAlleleExpression)
#library(Seurat)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(purrr)
library(stringi)
library(dplyr)
library(tidyr)



TestSNP_aod<-function(dat,SNPLevel=T,minCount=50,minSamp=10,form=cbind(alt, Num - alt) ~ 1,form2=cbind(alt, Num - alt) ~ 0,coef_ret="(Intercept)",numTest=-1)
{
    print("Format")
    if(SNPLevel)
    {
    dat<-dat %>% unite(Feature,Gene,SNP,sep="_",remove=F)
    }
    else{
        dat["Feature"]=dat[,"Gene"]
    }
    tab<-dat %>% group_by(Feature,Sample) %>% summarise(Num=length(unique(Allele)),Count=sum(Count)) %>% group_by(Feature) %>% summarise(Count=sum(Count),NumSamp=sum(Num>1)) %>% as.data.frame()
    print(length(unique(dat$Gene)))
    if(!SNPLevel)
    {
        minSamp=0;
    }
    feats=tab[tab$NumSamp>minSamp & tab$Count>minCount,"Feature"]
    dat=dat[dat$Feature %in% feats,]
    print("Number to Test:")
    print(length(feats))
    print(length(unique(dat$Gene)))



    print("Split by SNP")
    bySNP=split(dat,f=dat$Feature)
    if(numTest>0){bySNP=bySNP[sample(1:length(bySNP),numTest)]}
    print("Test")
    tic()
    nams=names(bySNP)
    out=lapply(names(bySNP),function(cur_feat){
        x=bySNP[[cur_feat]]
        x=x %>% spread(Allele,Count,fill=0)
        x["Num"]=x["alt"]+x["ref"]
        fit=tryCatch({betabin(form,~1,data=x)},error=function(cond){return(NULL)})
        if(is.null(fit))
        {
            return(NULL)
        }

        coef = tryCatch({summaryAOD(fit)@Coef},error = function(cond) {return(NULL)});
        if(is.null(coef))
        {
            return(NULL)
        }
        coef=data.frame(coef)
        colnames(coef)[4]="pval"
        coef["Test"]=rownames(coef)
        coef["SNP"]=cur_feat
        coef["NumSamp"]=length(unique(x$Sample))
        return(coef)

    })
    toc()
    bySNP=out
    names(bySNP)=nams
    bySNP[sapply(bySNP, is.null)]=NULL

    ret=do.call(rbind,bySNP)

    ret=data.frame(ret)


    ret=ret[order(ret$pval),]
    ret["logP"]=log(2)+pnorm(abs(ret[,"z.value"]),lower.tail=FALSE, log.p=TRUE)
    ret["pval"]=exp(ret[,"logP"])
    ret["padj"]=p.adjust(ret[,"pval"])
    ret["Gene"]=as.character(lapply(ret[,"SNP"],function(x){strsplit(x,split="_")[[1]][1]}))
    rownames(ret)=NULL
    return(ret)

}




dat=read.table("Z93S.QC.txt",header=T)
dat=dat[dat$total>10,]
for(i in c(2:11,15)){print(colnames(dat)[i]);dat[i]=dat[,i]/dat[,"total"]}
tab=read.table("/stanley/levin_asap_storage/612-eqtl/DS_ForPaper/GTExData/samp_Z93S_130/output/AlleleCounts/counts.txt")
seur=readRDS("../UMAP/seur.RDS")
seur=subset(seur,CellType=="Excitatory")
tab["Name"]=sub("^","Z93S_",tab[,1])
meta=seur@meta.data
tab=tab[tab$Name %in% rownames(meta),]
tab=tab[tab[,3]!="Ambig",]

colnames(tab)[2]="Gene"
tab["Sample"]="samp"
colnames(tab)[3]="Allele"
colnames(tab)[4]="Count"
tab["SNP"]=tab["Gene"]
tab[tab$Allele=="All1","Allele"]="alt"
tab[tab$Allele=="All2","Allele"]="ref"
num<-tab %>% group_by(Gene) %>% summarise(numAlt=sum(Allele=="alt"),numRef=sum(Allele=="ref")) %>% as.data.frame()
tab=tab[tab$Gene %in% num[num$numAlt>0 & num$numRef>0,"Gene"],]

rownames(dat)=dat[,1]
for(i in colnames(dat)){tab[i]=dat[tab[,1],i]}

out=map(colnames(dat)[2:16],function(i){
    print(i);tab["QC"]=as.numeric(scale(tab[,i]));if(sum(tab[,i])==0){print("Skip");return(NULL)};out=TestSNP_aod(tab,minCount=10,minSamp=0,form = cbind(alt,Num - alt) ~ 1 + QC,coef_ret = "QC")
})


names(out)=colnames(dat)[2:16]

saveRDS(out,"out.Z93S.RDS")




dat=read.table("YJ89.QC.txt",header=T)
dat=dat[dat$total>10,]
for(i in c(2:11,15)){print(colnames(dat)[i]);dat[i]=dat[,i]/dat[,"total"]}
tab=read.table("/stanley/levin_asap_storage/612-eqtl/DS_ForPaper/GTExData/samp_YJ89_130/output/AlleleCounts/counts.txt")
#seur=readRDS("../UMAP/seur.RDS")
#seur=subset(seur,CellType=="Excitatory")
tab["Name"]=sub("^","YJ89_",tab[,1])
meta=seur@meta.data
tab=tab[tab$Name %in% rownames(meta),]
tab=tab[tab[,3]!="Ambig",]

colnames(tab)[2]="Gene"
tab["Sample"]="samp"
colnames(tab)[3]="Allele"
colnames(tab)[4]="Count"
tab["SNP"]=tab["Gene"]
tab[tab$Allele=="All1","Allele"]="alt"
tab[tab$Allele=="All2","Allele"]="ref"
num<-tab %>% group_by(Gene) %>% summarise(numAlt=sum(Allele=="alt"),numRef=sum(Allele=="ref")) %>% as.data.frame()
tab=tab[tab$Gene %in% num[num$numAlt>0 & num$numRef>0,"Gene"],]

rownames(dat)=dat[,1]
for(i in colnames(dat)){tab[i]=dat[tab[,1],i]}

out=map(colnames(dat)[2:16],function(i){
    print(i);tab["QC"]=tab[,i];if(sum(tab[,i])==0){print("Skip");return(NULL)};out=TestSNP_aod(tab,minCount=10,minSamp=0,form = cbind(alt,Num - alt) ~ 1 + QC,coef_ret = "QC")
})

names(out)=colnames(dat)[2:16]

saveRDS(out,"out.YJ89.RDS")
