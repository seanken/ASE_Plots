library(scAlleleExpression)
#library(Seurat)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(purrr)
library(stringi)
library(dplyr)
library(tidyr)
library(qs)


peaks=read.table("/stanley/levin_asap_storage/612-eqtl/SingleNuc_data/UpdatedAIData/RunNewData/Sierra/MergeDS/merged.peaks.txt",header=T)

#meta=readRDS("/stanley/levin_asap_storage/612-eqtl/SingleNuc_data/UpdatedAIData/RunNewData/DataFromClemens/meta.with.files.and.Sierra.RDS")
meta=qread("/stanley/levin_asap_storage/ssimmons/ClemensUpdated/meta.with.Sierra.STARv2.8.qs")
tab=readRDS("../GTEx/celltype.interacting.eQTLs.SuppTab8.RDS")
head(meta)
meta=meta[meta$MajorCellTypes=="GLU_Neurons",]
tab=tab[tab$"Excitatory.BH-FDR"<.05,]
genes=tab[,2]

beta=tab[,"Excitatory.interaction.beta"]
tab["All2"]=map_chr(tab[,"SNP"],function(x)strsplit(x,"_")[[1]][2])
beta[tab$All2!=tab$Allele.assessed]=-beta[tab$All2!=tab$Allele.assessed]
snps_dirty=map_chr(tab[,3],function(x){s=strsplit(x,":")[[1]];chrom=s[1];pos=s[2];vals=s[4];vals=sub("_",":",vals);paste(chrom,pos,vals,sep=":")})
snps_dirty2=map_chr(tab[,3],function(x){s=strsplit(x,":")[[1]];chrom=s[1];pos=s[2];vals=s[4];vals=sub("_",":",vals);paste(chrom,pos,stri_reverse(vals),sep=":")})
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
tab_ie=tab

peaks=peaks[,c("Gene","polyA_ID")]
tab=inner_join(tab,peaks)

reload=T


snps=tab[,"SNP_name"]
genes=tab[,"polyA_ID"]
genes=gsub("-","_",gsub(":","_",genes,fixed=T),fixed=T)

#if(reload)
#{
dat=GetSNPs(meta=meta,genes=genes,snp_list=snps,samp="batchSamp",snp_col="SNP",cbc_col="CBC",cond="case",cellType="MajorCellTypes",bulk=T,all_col="Allele_Sierra")
saveRDS(dat,"counts.pseudo.RDS")
#}
#else{
#dat=readRDS("counts.pseudo.RDS")
#}

samps=unique(dat$Sample)
celltypes=unique(dat$CellType)
celltypes=celltypes[celltypes!="Unknown_Immune_Cells"]
out=map(rep((2:10)*9-1,1),function(x){
    print(x)
    set.seed(x)
    samps_rand=sample(samps,x)
    tab=dat[dat$Sample %in% samps_rand,]
    vals=map(celltypes,function(y){
        print(y)
        tmp=tab[tab$CellType==y,]
        mrk=tryCatch({TestSNP_aod(tmp,minCount = 10, minSamp = 5)},error=function(cond){print("Yuck!");return(NULL)})
        
        mrk["CellType"]=y
        return(mrk)
    })
    vals=vals %>% discard(is.null)
    vals=do.call(rbind,vals)
    vals["DS"]=x
    return(vals)
})




out1=do.call(rbind,out)
out=out1

saveRDS(out1,"AI.pseudo.RDS")

tab<-out %>% group_by(CellType,DS) %>% summarise(NumSig=sum(padj<.05)) %>% as.data.frame()
p=ggplot(tab,aes(x=DS,y=NumSig))+geom_point()+geom_line()+ylab("Number Signficant Peaks")+xlab("Number Individuals")
ggsave("DS.Sierra.pdf",p)


out=map(rep((2:10)*9-1,1),function(x){
    print(x)
    set.seed(x)
    samps_rand=sample(samps,x)
    tab=dat[dat$Sample %in% samps_rand,]
    vals=map(celltypes,function(y){
        print(y)
        tmp=tab[tab$CellType==y,]
        mrk=tryCatch({TestSNP_aod_Sierra(tmp,minCount = 10, minSamp = 5)},error=function(cond){print("Yuck!");return(NULL)})
        
        mrk["CellType"]=y
        return(mrk)
    })
    vals=vals %>% discard(is.null)
    vals=do.call(rbind,vals)
    vals["DS"]=x
    return(vals)
})




out2=do.call(rbind,out)
out=out2


saveRDS(out2,"AI.splicing.pseudo.RDS")

tab<-out %>% group_by(CellType,DS) %>% summarise(NumSig=sum(padj<.05,na.rm=T)) %>% as.data.frame()
p=ggplot(tab,aes(x=DS,y=NumSig))+geom_point()+geom_line()+ylab("Number Signficant Peaks (Splicing)")+xlab("Number Individuals")
ggsave("DS.Sierra.splice.pdf",p)

print("Compare to gene level")
out=readRDS("../GTEx/AI.pseudo.RDS")
out=out[out$DS==89,]
out1=out1[out1$DS==89,]
out["SNPname"]=map_chr(out[,"SNP"],function(x){s=strsplit(x,"_")[[1]];s[length(s)]})
out1["SNPname"]=map_chr(out1[,"SNP"],function(x){s=strsplit(x,"_")[[1]];s[length(s)]})
comb=inner_join(out,out1,by=c("SNPname","Gene"))
comb["Sig"]="Not Significant"
comb[comb$padj.x<.05,"Sig"]="Significant at gene level"
comb[comb$padj.y<.05,"Sig"]="Significant at peak level"
comb[comb$padj.y<.05 & comb$padj.x<.05,"Sig"]="Significant in both"
p=ggplot(comb,aes(x=Estimate.x,y=Estimate.y,color=Sig))+geom_point()+xlab("Effect size gene")+ylab("Effect size peak")+labs(color="Signficance")+geom_hline(yintercept=0,linetype="dotted")+geom_vline(xintercept=0,linetype="dotted")
ggsave("Gene.vs.peak.pdf",p,width=10)
