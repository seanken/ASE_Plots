library(scAlleleExpression)
#library(Seurat)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(purrr)

reload=FALSE

meta=readRDS("/stanley/levin_asap_storage/612-eqtl/SingleNuc_data/UpdatedAIData/RunNewData/DataFromClemens/meta.with.files.and.Sierra.RDS")
#dat=readRDS("/stanley/levin_asap_storage/612-eqtl/SingleNuc_data/UpdatedAIData/RunNewData/DataFromClemens/snps.allele.counts.table.full.RDS")
dat=readRDS("../GTEx/counts.pseudo.RDS")
#if(reload)
#{
#    genes=readRDS("/stanley/levin_asap_storage/612-eqtl/SingleNuc_data/UpdatedAIData/RunNewData/DataFromClemens/genes.all.RDS")
#    snps=readRDS("/stanley/levin_asap_storage/612-eqtl/SingleNuc_data/UpdatedAIData/RunNewData/DataFromClemens/snps.all.RDS")
#    dat=GetSNPs(meta=meta,genes=genes,snp_list=snps,samp="batchSamp",snp_col="SNP_Full",cbc_col="CBC",cond="case",cellType="MajorCellTypes",bulk=T)
#}

samps=unique(dat$Sample)
celltypes=unique(dat$CellType)
celltypes=celltypes[celltypes!="Unknown_Immune_Cells"]

out=map(1:100,function(x){
    print(x)
    tab=dat
    set.seed(x)
    
    for(samp in samps)
    {
        randVal=rbinom(1,1,.5)
        
        if(randVal==1)
        {
            
            tab[tab$Allele=="alt" & tab$Sample==samp,"Allele"]="ref2"
            tab[tab$Allele=="ref" & tab$Sample==samp,"Allele"]="alt"
            tab[tab$Allele=="ref2" & tab$Sample==samp,"Allele"]="ref"
        }
    }
    
    vals=map(celltypes,function(y){
        print(y)
        tmp=tab[tab$CellType==y,]
        mrk=tryCatch({TestSNP_aod(tmp,minCount = 10, minSamp = 5)},error=function(cond){print("Yuck!");return(NULL)})
        print(head(mrk))
        mrk["CellType"]=y
        return(mrk)
    })
    vals=vals %>% discard(is.null)
    vals=do.call(rbind,vals)
    vals["DS"]=x
    return(vals)
})




out1=do.call(rbind,out)
saveRDS(out1,"pseudo.results.RDS")


#dat=readRDS("counts.per.cell.RDS")
dat=readRDS("../GTEx/counts.cell.RDS")
#if(reload)
#{
#    genes=readRDS("/stanley/levin_asap_storage/612-eqtl/SingleNuc_data/UpdatedAIData/RunNewData/DataFromClemens/genes.all.RDS")
#    snps=readRDS("/stanley/levin_asap_storage/612-eqtl/SingleNuc_data/UpdatedAIData/RunNewData/DataFromClemens/snps.all.RDS")
#    dat=GetSNPs(meta=meta,genes=genes,snp_list=snps,samp="batchSamp",snp_col="SNP_Full",cbc_col="CBC",cond="case",cellType="MajorCellTypes",bulk=F)
#}
samps=unique(dat$Sample)
celltypes=unique(dat$CellType)
celltypes=celltypes[celltypes!="Unknown_Immune_Cells"]
out=map(1:100,function(x){
    print(x)
    tab=dat
    set.seed(x)
    for(samp in samps)
    {
        randVal=rbinom(1,1,.5)
        if(randVal==1)
        {
            tab[tab$Allele=="alt" & tab$Sample==samp,"Allele"]="ref2"
            tab[tab$Allele=="ref" & tab$Sample==samp,"Allele"]="alt"
            tab[tab$Allele=="ref2" & tab$Sample==samp,"Allele"]="ref"
        }
    }
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
saveRDS(out2,"mixedmodel.results.RDS")

out1["Type"]="Pseudobulk"
out2["Type"]="Mixed Model"

inter=intersect(colnames(out1),colnames(out2))

out=rbind(out1[,inter],out2[,inter])
out=out[!is.na(out$padj),]

tab<-out %>% group_by(Type,DS) %>% summarise(Perc=100*mean(pval<.05)) %>% as.data.frame()
ggplot(tab,aes(x=Type,y=Perc))+geom_boxplot()+geom_hline(yintercept=5,linetype="dotted")+ylab("Percent tests with p-value<.05")+xlab("Model used")
ggsave("FPR.pdf")

saveRDS(out,"DE.results.RDS")
out_gene=out

