library(DAESC)
library(dplyr);library(tidyr)

dat=readRDS("../GTEx/counts.cell.RDS")
dat <- dat %>% unite(Feature,Gene,SNP,remove=F,sep="_")
tab <- dat %>% group_by(Feature,Sample,Condition) %>% summarise(Num=sum(Count)) %>% group_by(Feature) %>% summarise(numSamp=length(Sample),Num=sum(Num)) %>% as.data.frame()
feats=tab[tab$numSamp>10 & tab$Num>1000,"Feature"]
dat=dat[dat$Feature %in% feats,]


 bySNP=split(dat,f=dat$Feature)

out=lapply(names(bySNP),function(cur_feat){
    print(cur_feat)
    x=bySNP[[cur_feat]]
    x<-x %>% spread(Allele,Count,fill=0)
    samps=unique(x$Sample[x$Condition=="PD"])
    x["Cond"]="ko"
    #numSamps=length(samps)
    #numSamps=floor(numSamps/2)
    #samps=sample(samps,numSamps)
    x[x$Sample %in% samps,"Cond"]="wt"
    x["n"]=x["ref"]+x["alt"]
    mod=model.matrix(~0+Cond,data=x)
    mod=as.matrix(mod[,"Condwt"])
    colnames(mod)="Condwt"
    res.bb <- tryCatch({daesc_mix(y=x$ref, n=x$n, subj=x$Sample, x=mod, niter=200, niter_laplace=2, num.nodes=3,optim.method="BFGS", converge_tol=1e-8)},error=function(cond){print("Yuck!");return(NULL)})
    print(cur_feat)
    print(res.bb)
    print(" ")
    return(res.bb)
})

saveRDS(out,"out.DAESC.mix.RDS")