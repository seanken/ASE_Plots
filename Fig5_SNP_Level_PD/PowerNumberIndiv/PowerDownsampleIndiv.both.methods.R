library(scAlleleExpression)
library(Seurat)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(purrr)
library(qs)

reload=FALSE


meta=qread('/stanley/levin_asap_storage/ssimmons/ClemensUpdated/meta.with.ASEInfo.qs')
dat=readRDS("../GTEx/counts.pseudo.RDS")
#meta=readRDS("/stanley/levin_asap_storage/612-eqtl/SingleNuc_data/UpdatedAIData/RunNewData/DataFromClemens/meta.with.files.and.Sierra.RDS")
#dat=readRDS("/stanley/levin_asap_storage/612-eqtl/SingleNuc_data/UpdatedAIData/RunNewData/DataFromClemens/snps.allele.counts.table.full.RDS")

if(reload)
{
    genes=readRDS("/stanley/levin_asap_storage/612-eqtl/SingleNuc_data/UpdatedAIData/RunNewData/DataFromClemens/genes.all.RDS")
    snps=readRDS("/stanley/levin_asap_storage/612-eqtl/SingleNuc_data/UpdatedAIData/RunNewData/DataFromClemens/snps.all.RDS")
    dat=GetSNPs(meta=meta,genes=genes,snp_list=snps,samp="batchSamp",snp_col="SNP_Full",cbc_col="CBC",cond="case",cellType="MajorCellTypes",bulk=T)
}

samps=unique(dat$Sample)
celltypes=unique(dat$CellType)
celltypes=celltypes[celltypes!="Unknown_Immune_Cells"]
out=map(rep((2:10)*9-1,1),function(x){
    print(x)
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



dat=readRDS("../GTEx/counts.per.cell.RDS")
if(reload)
{
    genes=readRDS("/stanley/levin_asap_storage/612-eqtl/SingleNuc_data/UpdatedAIData/RunNewData/DataFromClemens/genes.all.RDS")
    snps=readRDS("/stanley/levin_asap_storage/612-eqtl/SingleNuc_data/UpdatedAIData/RunNewData/DataFromClemens/snps.all.RDS")
    dat=GetSNPs(meta=meta,genes=genes,snp_list=snps,samp="batchSamp",snp_col="SNP_Full",cbc_col="CBC",cond="case",cellType="MajorCellTypes",bulk=F)
}
samps=unique(dat$Sample)
celltypes=unique(dat$CellType)
celltypes=celltypes[celltypes!="Unknown_Immune_Cells"]
out=map(rev(rep((2:10)*9-1,1)),function(x){
    print(x)
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

out=rbind(out1,out2)

tab<-out %>% group_by(CellType,DS,Gene,Type) %>% summarise(NumSig=sum(pval<.05/89),Num=length(Gene)) %>% group_by(CellType,DS,Type) %>% summarise(Len=max(Num),NumSig=sum(NumSig)/Len) %>% as.data.frame()
p=ggplot(tab,aes(x=DS,y=NumSig,color=CellType))+geom_point()+geom_line()+xlab("Average Number of Individuals Used")+ylab("Number of Signficant Genes")+facet_wrap(~Type)
ggsave("NumberSig.vs.NumberIndiv.pdf")
saveRDS(out,"DE.results.RDS")
out_gene=out

print("Run For Sierra")
dat=readRDS("....")
if(reload)
{
    genes=readRDS("/stanley/levin_asap_storage/612-eqtl/SingleNuc_data/UpdatedAIData/RunNewData/DataFromClemens/genes.all.RDS")
    snps=readRDS("/stanley/levin_asap_storage/612-eqtl/SingleNuc_data/UpdatedAIData/RunNewData/DataFromClemens/snps.all.RDS")
    dat=GetSNPs(meta=meta,genes=genes,snp_list=snps,samp="batchSamp",snp_col="SNP_Full",cbc_col="CBC",cond="case",cellType="MajorCellTypes",bulk=T,all_col="Allele_Sierra")
}
samps=unique(dat$Sample)
celltypes=unique(dat$CellType)
celltypes=celltypes[celltypes!="Unknown_Immune_Cells"]
out=map(rev(rep((2:10)*9-1,1)),function(x){
    print(x)
    samps_rand=sample(samps,x)
    tab=dat[dat$Sample %in% samps_rand,]
    vals=map(celltypes,function(y){
        print(y)
        tmp=tab[tab$CellType==y,]
        mrk=tryCatch({TestSNP_aod(tmp,minCount = 10, minSamp = 5)},error=function(cond){print("Yuck!");return(NULL)})
        print(sum(mrk$pval<.05/89))
        mrk["CellType"]=y
        return(mrk)
    })
    vals=vals %>% discard(is.null)
    vals=do.call(rbind,vals)
    vals["DS"]=x
    return(vals)
})

out=do.call(rbind,out)
saveRDS(out,"DE.Sierra.results.RDS")