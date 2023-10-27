library(scAlleleExpression)
#library(Seurat)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(purrr)
library(stringi)
library(dplyr)
library(tidyr)

reload=FALSE
#Run "cat /stanley/levin_asap_storage/612-eqtl/SingleNuc_data/UpdatedAIData/RunNewData/samp_*/output/SNPLevelCounts/comb.bed | awk '{print $1"\t"$2}' | sort | uniq > snp.to.genes.txt" first
print("Load")
args = commandArgs(trailingOnly=TRUE)

meta=readRDS("/stanley/levin_asap_storage/612-eqtl/SingleNuc_data/UpdatedAIData/RunNewData/DataFromClemens/meta.with.files.and.Sierra.RDS")

tab=readRDS("forTrans.RDS")
tab_all=readRDS("ASE.for.trans.RDS")
#genes=tab[,"Gene.symbol"]
if(reload)
{
    genes=tab[,"Gene.symbol"]
    tab=readRDS("trans.100Pcs.nocrossmap.RDS")
    beta=map_dbl(tab[,"Meta-analysis.Beta.and.SE"],function(x) as.numeric(strsplit(x," ")[[1]][1]))
    tab["All2"]=map_chr(tab[,"SNP"],function(x)strsplit(x,"_")[[1]][2])
    beta[tab$All2!=tab$Effect.Allele]=-beta[tab$All2!=tab$Effect.Allele]
    snps_dirty=map_chr(tab[,2],function(x){s=strsplit(x,":")[[1]];chrom=s[1];pos=s[2];vals=s[4];vals=sub("_",":",vals);paste(chrom,pos,vals,sep=":")})
    snps_dirty2=map_chr(tab[2],function(x){s=strsplit(x,":")[[1]];chrom=s[1];pos=s[2];vals=s[4];vals=sub("_",":",vals);paste(chrom,pos,stri_reverse(vals),sep=":")})
    snps_all=system("zcat /stanley/levin_asap_storage/612-eqtl/GenotypeData/Clean_vcf/Combined/comb_new.no.chr.vcf.gz | awk '{print $3}'",intern=T)
    vals=snps_dirty %in% snps_all
    snps=map_chr(1:length(snps_dirty),function(x){if(vals[x]){snps_dirty[x];}else{snps_dirty2[x]}})
    beta=map_dbl(1:length(snps_dirty),function(x){if(vals[x]){beta[x];}else{-beta[x]}})
    vals=snps %in% snps_all
    snps=snps[vals]
    genes=genes[vals]
    beta=beta[vals]
    tab=data.frame(SNP_name=snps,Gene=genes,beta=beta)
    tab<-tab %>% unite(SNP,Gene,SNP_name,remove=F,sep="_")


    print("Get cis Gene")
    snp2gene=read.table("snp.to.genes.txt")
    colnames(snp2gene)=c("SNP_name","Gene_cis")
    tab=inner_join(tab,snp2gene)
    #out=map(meta[,"SNP"],function(x){tmp=read.table(x);})
    #out=do.call(rbind,out)
    #colnames(out)=c("SNP","Gene","allele")
    #out<-out %>% group_by(SNP,Gene) %>% summarise() %>% as.data.frame()
    #...
    saveRDS(tab,"forTrans.RDS")

    print("Load Expression")

    print("Load data")
    genes=tab[,"Gene_cis"]
    snps=tab[,"SNP_name"]
    tab_all=GetSNPs(meta=meta,genes=genes,snp_list=snps,samp="batchSamp",snp_col="SNP_Full",cbc_col="CBC",cond="case",cellType="MajorCellTypes",bulk=F)
    saveRDS(tab_all,"ASE.for.trans.RDS")
}
dat=readRDS("/stanley/levin_asap_storage/612-eqtl/SingleNuc_data/UpdatedAIData/RunNewData/DataFromClemens/TransEQTL/expression.1perc.RDS")

print("Test!")
#res=apply(tab_all,1,function(x){
pos=args[1]
#x=as.character(tab_all[pos,])
gene=as.character(tab[pos,"Gene_cis"]);
snp=as.character(tab[pos,"SNP_name"]);
print(gene)
print(snp)
print("")
#GetTransEQTL(meta,gene,snp,tab_all,dat,snp_col="SNP_Full",samp_col="batchSamp",raw_col="dir_ASE")
res=tryCatch({GetTransEQTL(meta,gene,snp,tab_all,dat,snp_col="SNP_Full",samp_col="batchSamp",raw_col="dir_ASE")},error=function(x){print("Yuck");return(NULL)})
#})
system("mkdir output")
saveRDS(res,paste("output/Trans.results.",gene,".",gsub(":","_",snp),".RDS",sep=""))