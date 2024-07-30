library(scAlleleExpression)
library(purrr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(DAESC)
library(qs)
source("RunDAESC.R")

dat_pb=qread("../Selection/dat1.qs") 
#meta=qread("meta.select.qs")  
dat_sn=qread("../Selection/scdat1.qs")
dat_pb=dat_pb[dat_pb$Condition!="ILB",]
dat_sn=dat_sn[dat_sn$Condition!="ILB",]

celltypes=unique(dat_sn$CellType)[1:8]
#celltypes=celltypes[celltypes!="Unknown_Immune_Cells"]


print("pseudobulk")
vals_pb=map(celltypes,function(y){
    print(y)
    tmp=dat_pb[dat_pb$CellType==y,]
    tic()
    mrk=tryCatch({TestSNP_aod(tmp,minCount = 10, minSamp = 5,form=cbind(alt, Num - alt) ~ Condition)},error=function(cond){print("Yuck!");return(NULL)})
    ret=toc()    
    mrk["CellType"]=y
    mrk=mrk[mrk$Test!="(Intercept)",]
    mrk["padj"]=p.adjust(mrk[,"pval"],"fdr")
    ret[["AI"]]=mrk
    return(ret)
})
names(vals_pb)=celltypes
qsave(vals_pb,"vals_pb.qs")


print("mixed model")
vals_sn=map(celltypes,function(y){
    print(y)
    tmp=dat_sn[dat_sn$CellType==y,]
    tic()
    mrk=tryCatch({TestSNP_mixed(tmp,minCount = 10, minSamp = 5,form=cbind(alt, Num - alt) ~ Condition + (1| Sample))},error=function(cond){print("Yuck!");return(NULL)})
    ret=toc()    
    mrk["CellType"]=y
    mrk=mrk[mrk$Test!="(Intercept)",]
    mrk["padj"]=p.adjust(mrk[,"pval"],"fdr")
    ret[["AI"]]=mrk
    return(ret)
})
names(vals_sn)=celltypes
qsave(vals_sn,"vals_sn.qs")


print("DAESC")
vals_ds=map(celltypes,function(y){
    print(y)
    tmp=dat_sn[dat_sn$CellType==y,]
    tic()
    mrk=tryCatch({mrk=TestSNP_DAESC(tmp,minCount = 10, minSamp = 5,form= ~ Condition)},error=function(cond){print("Yuck outer!");print(cond);return(NULL)})
    ret=toc()    
    mrk["CellType"]=y
    #mrk=mrk[mrk$Test!="(Intercept)",]
    #if(!is.null())mrk["padj"]=p.adjust(mrk[,"pval"],"fdr")
    ret[["AI"]]=mrk
    return(ret)
})
names(vals_ds)=celltypes
qsave(vals_ds,"vals_ds.qs")
