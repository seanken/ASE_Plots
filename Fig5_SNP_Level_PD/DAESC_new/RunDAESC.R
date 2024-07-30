TestSNP_DAESC<-function(dat,minCount = 10, minSamp = 5,form= ~ Condition,form2=~1,numTest=0,coef=2)
{
    print("Format")
    dat<-dat %>% unite(Feature,Gene,SNP,sep="_",remove=F)
    tab<-dat %>% group_by(Feature,Sample) %>% summarise(Num=length(unique(Allele)),Count=sum(Count)) %>% group_by(Feature) %>% summarise(Count=sum(Count),NumSamp=sum(Num>1)) %>% as.data.frame()
    print(length(unique(dat$Gene)))
    feats=tab[tab$NumSamp>minSamp & tab$Count>minCount,"Feature"]
    dat=dat[dat$Feature %in% feats,]
    print("Number to Test:")
    print(length(feats))
    print(length(unique(dat$Gene)))
    bySNP=split(dat,f=dat$Feature)
    if(numTest>0){bySNP=bySNP[1:numTest]}
    print("Test")
    tic()
    nams=names(bySNP)

    out=lapply(names(bySNP),function(cur_feat){
        #print(cur_feat)
        x=bySNP[[cur_feat]]
        x<-x %>% spread(Allele,Count,fill=0)
        x["n"]=x["ref"]+x["alt"]
        mod=model.matrix(~1+Condition,data=x)
        #mod=as.matrix(mod[,"ConditionPD"])
        #colnames(mod)="ConditionPD"
        #print(head(mod))
        #print(head(x))
        res.bb <- tryCatch({daesc_bb(y=x$ref, n=x$n, subj=x$Sample, x=mod, niter=200, niter_laplace=2, num.nodes=3,optim.method="BFGS", converge_tol=1e-8)},error=function(cond){return(NULL)})
        #print(cur_feat)
        #print(res.bb)
        #print(" ")
        return(res.bb)
    })

    print(out)
    bySNP=out
    names(bySNP)=nams
    bySNP[sapply(bySNP, is.null)]=NULL
    mrk=data.frame(Test=names(bySNP),Estimate=map_dbl(bySNP,function(x) x$b[coef]),pval=map_dbl(bySNP,function(x) x$p.value),Note=map_chr(bySNP,function(x) x$note),NoteNull=map_chr(bySNP,function(x) x$note.null))
    print(head(mrk))
    mrk["padj"]=p.adjust(mrk[,"pval"],"fdr")
    mrk=mrk[order(mrk$pval),]
    return(mrk)

}
