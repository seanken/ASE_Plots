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

args = commandArgs(trailingOnly=TRUE)
num=as.numeric(args[1])

dat_pb=qread("../Selection/dat1.qs") 
#meta=qread("meta.select.qs")  
dat_sn=qread("../Selection/scdat1.qs")
dat_pb=dat_pb[dat_pb$Condition!="ILB",]
dat_sn=dat_sn[dat_sn$Condition!="ILB",]

#celltypes=unique(dat_pb$CellType)
#celltypes=celltypes[celltypes!="Unknown_Immune_Cells"]
celltypes=unique(dat_sn$CellType)[1:8]



print("DAESC")
print(num)
y=celltypes[num]
print(celltypes)
#vals_ds=map(celltypes,function(y){
    print(y)
    tmp=dat_sn[dat_sn$CellType==y,]
    tic()
    mrk=tryCatch({mrk=TestSNP_DAESC(tmp,minCount = 10, minSamp = 5,form= ~ Condition)},error=function(cond){print("Yuck outer!");print(cond);return(NULL)})
    ret=toc()    
    mrk["CellType"]=y
    #mrk=mrk[mrk$Test!="(Intercept)",]
    #if(!is.null())mrk["padj"]=p.adjust(mrk[,"pval"],"fdr")
    ret[["AI"]]=mrk
#    return(ret)
#})
#ret
#names(vals_ds)=celltypes
qsave(ret,paste0("vals_ds.",y,".qs"))
