library(DAESC)
library(dplyr);library(tidyr)

out=readRDS("/stanley/levin_asap_storage/612-eqtl/SingleNuc_data/UpdatedAIData/RunNewData/DataFromClemens/gene.level.split.RDS")
args = commandArgs(trailingOnly=TRUE)
num=args[1]
out=out[[num]]
dat=out

bySNP=split(dat,f=dat$Gene)

out=lapply(names(bySNP),function(cur_feat){
    print(cur_feat)
    x=bySNP[[cur_feat]]
    x["Condition"]=x[,"Cond"]
    x=x[x$CellType=="GLU_Neurons",]
    if(dim(x)[1]<50)
    {
        return(NULL)
    }
    
    x<-x %>% spread(Allele,Count,fill=0)
    samps=unique(x$Sample[x$Condition=="PD"])
    #x["Cond"]="ko"
    #numSamps=length(samps)
    #numSamps=floor(numSamps/2)
    #samps=sample(samps,numSamps)
    #x[x$Sample %in% samps,"Cond"]="wt"
    x["n"]=x["ref"]+x["alt"]
    mod1=model.matrix(~1+Condition,data=x)
    mod2=model.matrix(~1,data=x)
    #mod=as.matrix(mod[,"Condwt"])
    print(head(mod1))
    print(head(mod2))
    #colnames(mod)="Condwt"
    res.bb <- tryCatch({daesc_mix(y=x$ref, n=x$n, subj=x$Sample, x=mod1,xnull=mod2, niter=200, niter_laplace=2, num.nodes=3,optim.method="BFGS", converge_tol=1e-8)},error=function(cond){print("Yuck!");return(NULL)})
    print(cur_feat)
    print(res.bb)
    print(" ")
    return(res.bb)
})
names(out)=names(bySNP)
system("mkdir output")
saveRDS(out,paste("output/out.",num,".RDS",sep=""))
