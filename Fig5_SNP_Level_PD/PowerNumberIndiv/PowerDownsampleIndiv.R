library(scAlleleExpression)
library(Seurat)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(purrr)

meta=readRDS("/stanley/levin_asap_storage/612-eqtl/SingleNuc_data/UpdatedAIData/RunNewData/DataFromClemens/meta.with.files.and.Sierra.RDS")
dat=readRDS("/stanley/levin_asap_storage/612-eqtl/SingleNuc_data/UpdatedAIData/RunNewData/DataFromClemens/snps.allele.counts.table.full.RDS")
samps=unique(dat$Sample)
celltypes=unique(dat$CellType)
celltypes=celltypes[celltypes!="Unknown_Immune_Cells"]
out=map(rep((2:10)*9-1,50),function(x){
    print(x)
    tab=dat[dat$Sample %in% sample(samps,x),]
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

out=do.call(rbind,out)
tab<-out %>% group_by(CellType,DS,Gene) %>% summarise(NumSig=sum(pval<.05/89),Num=length(Gene)) %>% group_by(CellType,DS) %>% summarise(Len=max(Num),NumSig=sum(NumSig)/Len) %>% as.data.frame()
p=ggplot(tab,aes(x=DS,y=NumSig,color=CellType))+geom_point()+geom_line()+xlab("Average Number of Individuals Used")+ylab("Number of Signficant Genes")
ggsave("NumberSig.vs.NumberIndiv.pdf")