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


reload=T

meta=qread('/stanley/levin_asap_storage/ssimmons/ClemensUpdated/meta.with.ASEInfo.qs')
#meta=readRDS("/stanley/levin_asap_storage/612-eqtl/SingleNuc_data/UpdatedAIData/RunNewData/DataFromClemens/meta.with.files.and.Sierra.RDS")
tab=readRDS("celltype.interacting.eQTLs.SuppTab8.RDS")
head(meta)
meta=meta[meta$MajorCellTypes=="GLU_Neurons",]
# tab=tab[tab$"Excitatory.BH-FDR"<.05,]
# genes=tab[,2]

# beta=tab[,"Excitatory.interaction.beta"]
# tab["All2"]=map_chr(tab[,"SNP"],function(x)strsplit(x,"_")[[1]][2])
# beta[tab$All2!=tab$Allele.assessed]=-beta[tab$All2!=tab$Allele.assessed]
# snps_dirty=map_chr(tab[,3],function(x){s=strsplit(x,":")[[1]];chrom=s[1];pos=s[2];vals=s[4];vals=sub("_",":",vals);paste(chrom,pos,vals,sep=":")})
# snps_dirty2=map_chr(tab[,3],function(x){s=strsplit(x,":")[[1]];chrom=s[1];pos=s[2];vals=s[4];vals=sub("_",":",vals);paste(chrom,pos,stri_reverse(vals),sep=":")})
# snps_all=system("zcat /stanley/levin_asap_storage/612-eqtl/GenotypeData/Clean_vcf/Combined/comb_new.no.chr.vcf.gz | awk '{print $3}'",intern=T)
# vals=snps_dirty %in% snps_all
# snps=map_chr(1:length(snps_dirty),function(x){if(vals[x]){snps_dirty[x];}else{snps_dirty2[x]}})
# beta=map_dbl(1:length(snps_dirty),function(x){if(vals[x]){beta[x];}else{-beta[x]}})
# vals=snps %in% snps_all
# snps=snps[vals]
# genes=genes[vals]
# beta=beta[vals]
# tab=data.frame(SNP_name=snps,Gene=genes,beta=beta)
# tab<-tab %>% unite(SNP,Gene,SNP_name,remove=F,sep="_")
# tab_ie=tab



dat=NULL




#if(reload)
#{
#dat=GetSNPs(meta=meta,genes=genes,snp_list=snps,samp="batchSamp",snp_col="SNP",cbc_col="CBC",cond="case",cellType="MajorCellTypes",bulk=T)
#saveRDS(dat,"counts.pseudo.RDS")
#}
#else{
dat=readRDS("counts.pseudo.RDS")
#}

samps=unique(dat$Sample)
celltypes=unique(dat$CellType)
celltypes=celltypes[celltypes!="Unknown_Immune_Cells"]
print("DE Start!")
out=map(c(91,rep((2:10)*9-1,1)),function(x){
    print("New Rep")
    print(x)
    set.seed(x)
    samps_rand=sample(samps,x)
    tab=dat[dat$Sample %in% samps_rand,]
    vals=map(celltypes,function(y){
        print("Cell Type")
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


print("DE Done!")

out1=do.call(rbind,out)
out=out1
comb=inner_join(out[out$DS==89,],tab_ie)
p=ggplot(comb[order(comb$padj,decreasing=T),],aes(x=beta,y=Estimate,color=padj<.05))+geom_point()+scale_color_manual(values=c("grey","red"))+geom_hline(yintercept=0,linetype="dotted")+geom_vline(xintercept=0,linetype="dotted")+ylab("Effect Size Single Nucleus Allelic Imbalance")+xlab("Bulk eQTL exctiatory neuron interaction coefficient")
ggsave("Corr.bulk.interaction.vs.ai.pdf",p,width=14)

comb=inner_join(out[out$DS==89,],tab_ie)
p=ggplot(comb[order(comb$padj,decreasing=T),],aes(x=beta,y=Estimate,color=padj<.05))+geom_point()+scale_color_manual(values=c("grey","red"))+geom_hline(yintercept=0,linetype="dotted")+geom_vline(xintercept=0,linetype="dotted")+ylab("Bulk eQTL Meta Coefficient")
ggsave("Corr.bulk.meta.vs.ai.pdf",p,width=14)


saveRDS(out1,"AI.pseudo.RDS")
print("Start mixed")
#dat=NULL
#if(reload)
#{
#dat=GetSNPs(meta=meta,genes=genes,snp_list=snps,samp="batchSamp",snp_col="SNP",cbc_col="CBC",cond="case",cellType="MajorCellTypes",bulk=F)
#saveRDS(dat,"counts.cell.RDS")
#}
#else{
dat=readRDS("counts.cell.RDS")
#}
samps=unique(dat$Sample)
celltypes=unique(dat$CellType)
celltypes=celltypes[celltypes!="Unknown_Immune_Cells"]
print(celltypes)
out=map(c(91,rev(rep((2:10)*9-1,1))),function(x){
    print(x)
    set.seed(x)
    samps_rand=sample(samps,x)
    tab=dat[dat$Sample %in% samps_rand,]
    vals=map(celltypes,function(y){
        print(y)
        tmp=tab[tab$CellType==y,]
        mrk=tryCatch({TestSNP_mixed(tmp,minCount = 10, minSamp = 5)},error=function(cond){print("Yuck!");return(NULL)})
        print(sum(mrk$pval<.05/89))
        mrk["CellType"]=y
        return(mrk)
    })
    vals=vals %>% discard(is.null)
    vals=do.call(rbind,vals)
    vals["DS"]=x
    return(vals)
})

out2=do.call(rbind,out)

out1["Type"]="Pseudobulk"
out2["Type"]="Mixed Model"
saveRDS(out2,"AI.mixed.RDS")
inter=intersect(colnames(out1),colnames(out2))
out=rbind(out1[,inter],out2[,inter])
out=out[!is.na(out$padj),]
saveRDS(out,"results.test.RDS")

tab<-out %>% group_by(CellType,DS,Gene,Type) %>% summarise(NumSig=sum(pval<.05/89),Num=length(Gene)) %>% group_by(CellType,DS,Type) %>% summarise(Len=max(Num),NumSig=sum(NumSig)/Len) %>% as.data.frame()
p=ggplot(tab,aes(x=DS,y=NumSig,color=Type))+geom_point()+geom_line()+xlab("Average Number of Individuals Used")+ylab("Number of Signficant Genes") #+facet_wrap(~Type)
ggsave("NumberSig.vs.NumberIndiv.pdf",p)
saveRDS(out,"DE.results.RDS")

#out=out[out$DS==89,]
comb=inner_join(out[out$DS==91,],tab_ie)
p=ggplot(comb[order(comb$padj,decreasing=T),],aes(x=beta,y=Estimate,color=padj<.05))+geom_point()+scale_color_manual(values=c("grey","red"))+geom_hline(yintercept=0,linetype="dotted")+geom_vline(xintercept=0,linetype="dotted")+ylab("Effect Size Single Nucleus Allelic Imbalance")+xlab("Bulk eQTL exctiatory neuron interaction coefficient")+facet_wrap(~Type)
ggsave("Corr.bulk.interaction.vs.ai.both.pdf",p,width=14,height=5)
