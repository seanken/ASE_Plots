library(scAlleleExpression)
library(purrr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(qs)

print("Set up")

#meta1=qread("meta.select.DS.qs")
meta1=qread("/stanley/levin_asap_storage/ssimmons/ClemensUpdated/meta.selection.without.rep.cells.DS.qs")
meta2=qread("meta.no.select.qs")
genes=qread("genes.qs")
snps=qread("snps.qs")

inter=intersect(rownames(meta1),rownames(meta2))
meta1=meta1[inter,]
meta2=meta2[inter,]

print("Coverage plot")

tab1<-meta1 %>% group_by(batchSamp,dirASE) %>% summarise(Num=length(age)) %>% as.data.frame()
tab1["nReads"]=map_dbl(tab1[,"dirASE"],function(x) read.csv(paste0(x,"/output/STARSolo/output/resultsSolo.out/GeneFull/Summary.csv"),header=F)[1,2])
tab2<-meta2 %>% group_by(batchSamp,dirASE) %>% summarise(Num=length(age)) %>% as.data.frame()
tab2["nReads"]=map_dbl(tab2[,"dirASE"],function(x) read.csv(paste0(x,"/output/STARSolo/output/resultsSolo.out/GeneFull/Summary.csv"),header=F)[1,2])
comb=inner_join(tab1,tab2,by="batchSamp")


p=ggplot(comb,aes(x=nReads.x/Num.x,y=nReads.y/Num.x))+geom_point()+geom_abline(slope=1,linetype="dotted")+ylim(c(0,160000))+ylab("Reads per cell, no selection")+xlab("Read per cell, selection")
ggsave("Reads.Per.Cell.Per.Sample.DS.pdf")
qsave(comb,"comb.DS.qs")

print("Load ASE data")
print("Selected")
dat1=GetSNPs(meta=meta1,genes=genes,snp_list=snps,samp="batchSamp",snp_col="SNP",cbc_col="CBC",cond="case",cellType="MajorCellTypes",bulk=T)
qsave(dat1,"dat1.DS.qs")
print("Not Selected")
dat2=GetSNPs(meta=meta2,genes=genes,snp_list=snps,samp="batchSamp",snp_col="SNP",cbc_col="CBC",cond="case",cellType="MajorCellTypes",bulk=T)
qsave(dat2,"dat2.DS.qs")


print("Test")
samps=unique(dat1$Sample)
celltypes=unique(dat1$CellType)
celltypes=celltypes[celltypes!="Unknown_Immune_Cells"]

vals1=map(celltypes,function(y){
    print(y)
    tab=dat1
    tmp=tab[tab$CellType==y,]
    mrk=tryCatch({TestSNP_aod(tmp,minCount = 10, minSamp = 5)},error=function(cond){print("Yuck!");return(NULL)})
        
    mrk["CellType"]=y
    return(mrk)
})
qsave(vals1,"vals1.DS.qs")
vals2=map(celltypes,function(y){
    print(y)
    tab=dat2
    tmp=tab[tab$CellType==y,]
    mrk=tryCatch({TestSNP_aod(tmp,minCount = 10, minSamp = 5)},error=function(cond){print("Yuck!");return(NULL)})
        
    mrk["CellType"]=y
    return(mrk)
})

qsave(vals2,"vals2.DS.qs")
