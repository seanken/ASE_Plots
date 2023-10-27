library(scAlleleExpression)
library(Seurat)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(purrr)



print("Load Seurat Data")
seur=readRDS("../UMAP/seur.RDS")

print("Set up data frame")
counts=system("ls count/*txt",intern=T)
ASE=map(counts,function(x){
    tab=read.table(x)
    samp=sub("count/samp_","",sub(".count.txt","",x))
    tab["Name"]=sub("^",paste(samp,"_",sep=""),tab[,1])
    print(dim(tab))
    tab=tab[tab$Name %in% names(seur@active.ident),]
    print(dim(tab))
    colnames(tab)[2]="Gene"
    tab=tab[grep(",",tab$Gene,invert=T),]
    tab["Sample"]="samp"
    colnames(tab)[3]="Allele"
    colnames(tab)[4]="Count"
    tab["SNP"]=tab["Gene"]
    tab[tab$Allele=="All1","Allele"]="alt"
    tab[tab$Allele=="All2","Allele"]="ref"
    num<-tab %>% group_by(Gene) %>% summarise(numAlt=sum(Allele=="alt"),numRef=sum(Allele=="ref")) %>% as.data.frame()
    tab=tab[tab$Gene %in% num[num$numAlt>0 & num$numRef>0,"Gene"],]
    out=TestSNP_aod(tab,minCount=10,minSamp=0)
    return(out)
})

saveRDS(ASE,"ASE.RDS")


TestSNP_aod_Sierra<-function (dat, minCount = 50, minSamp = 10, form = cbind(Count, 
    Tot - Count) ~ Allele + Sample, form2 = cbind(alt, Num - 
    alt) ~ 0, coef_ret = "Allele", numTest = -1) 
{
    print("Format")
    dat <- dat %>% unite(Feature, Gene, SNP, sep = "_", remove = F)
    tab <- dat %>% group_by(Feature, Sample) %>% summarise(Num = length(unique(Allele)), 
        Count = sum(Count)) %>% group_by(Feature) %>% summarise(Count = sum(Count), 
        NumSamp = sum(Num > 1)) %>% as.data.frame()
    print(length(unique(dat$Gene)))
    feats = tab[tab$NumSamp > minSamp & tab$Count > minCount, 
        "Feature"]
    dat = dat[dat$Feature %in% feats, ]
    print("Number to Test:")
    print(length(feats))
    print(length(unique(dat$Gene)))
    dat["GeneName"] = as.character(lapply(dat[, "Gene"], function(x) {
        strsplit(x, "_")[[1]][1]
    }))
    tab <- dat %>% group_by(Allele, GeneName, Sample) %>% summarise(Tot = sum(Count)) %>% 
        as.data.frame()
    dat <- dat %>% spread(Allele, Count, fill = 0) %>% gather(Allele, 
        Count, alt, ref)
    dat = inner_join(dat, tab)
    print(head(dat))
    print("Split by SNP")
    bySNP = split(dat, f = dat$Feature)
    if (numTest > 0) {
        bySNP = bySNP[sample(1:length(bySNP), numTest)]
    }
    print("Test")
    tic()
    nams = names(bySNP)
    out = lapply(names(bySNP), function(cur_feat) {
        x = bySNP[[cur_feat]]
        fit = tryCatch({
            betabin(form, ~1, data = x)
        }, error = function(cond) {
            return(NULL)
        })
        if (is.null(fit)) {
            return(NULL)
        }
        coef = tryCatch({
            summaryAOD(fit)@Coef
        }, error = function(cond) {
            return(NULL)
        })
        if (is.null(coef)) {
            return(NULL)
        }
        coef = data.frame(coef)
        colnames(coef)[4] = "pval"
        coef["Test"] = rownames(coef)
        coef = coef[grep(coef_ret, rownames(coef)), ]
        coef["SNP"] = cur_feat
        coef["NumSamp"] = length(unique(x$Sample))
        return(coef)
    })
    toc()
    bySNP = out
    names(bySNP) = nams
    bySNP[sapply(bySNP, is.null)] = NULL
    ret = do.call(rbind, bySNP)
    ret = data.frame(ret)
    ret = ret[order(ret$pval), ]
    ret["logP"] = log(2) + pnorm(abs(ret[, "z.value"]), lower.tail = FALSE, 
        log.p = TRUE)
    ret["pval"] = exp(ret[, "logP"])
    ret["padj"] = p.adjust(ret[, "pval"], "fdr")
    ret["Gene"] = as.character(lapply(ret[, "SNP"], function(x) {
        strsplit(x, split = "_")[[1]][1]
    }))
    rownames(ret) = NULL
    return(ret)
}



ASE2=map(counts,function(x){
    tab=read.table(x)
    samp=sub("count/samp_","",sub(".count.txt","",x))
    tab["Name"]=sub("^",paste(samp,"_",sep=""),tab[,1])
    print(dim(tab))
    tab=tab[tab$Name %in% names(seur@active.ident),]
    print(dim(tab))
    colnames(tab)[2]="Gene"
    tab=tab[grep(",",tab$Gene,invert=T),]
    tab["Sample"]="samp"
    colnames(tab)[3]="Allele"
    colnames(tab)[4]="Count"
    tab["SNP"]=tab["Gene"]
    tab[tab$Allele=="All1","Allele"]="alt"
    tab[tab$Allele=="All2","Allele"]="ref"
    num<-tab %>% group_by(Gene) %>% summarise(numAlt=sum(Allele=="alt"),numRef=sum(Allele=="ref")) %>% as.data.frame()
    tab=tab[tab$Gene %in% num[num$numAlt>0 & num$numRef>0,"Gene"],]
    out=TestSNP_aod_Sierra(tab,minCount=10,minSamp=0,form = cbind(Count,Tot - Count) ~ Allele)
    return(out)
})

saveRDS(ASE2,"ASE2.RDS")

genes1=map(ASE,function(x) unique(x[x$padj<.05 & !is.na(x$padj),"Gene"]))
genes2=map(ASE2,function(x) unique(x[x$padj<.05 & !is.na(x$padj),"Gene"]))


lst=system("ls /stanley/levin_asap_storage/612-eqtl/DS_ForPaper/GTExData/samp_* | grep -v Hyb | grep : | sed 's/://g' | grep '0$\\|5$' | grep 130 ",intern=T)
counts=sub("$","/output/AlleleCounts/counts.txt",lst)

len=map_dbl(lst,function(x){s=strsplit(x,"_")[[1]];as.numeric(s[length(s)])})
samp=map_chr(lst,function(x){s=strsplit(x,"_")[[1]];s[length(s)-1]})

fils=data.frame("fil"=counts,"len"=len,"samp"=samp)
fils["QC"]=counts=sub("$","/output/STARSolo/output/resultsSolo.out/GeneFull/Summary.csv",lst)


ASE_gene=apply(fils,1,function(x){
    print(x["len"])
    tab=read.table(x["fil"])
    tab["Name"]=sub("^",paste(x["samp"],"_",sep=""),tab[,1])
    print(dim(tab))
    tab=tab[tab$Name %in% names(seur@active.ident),]
    print(dim(tab))
    colnames(tab)[2]="Gene"
    tab["Sample"]="samp"
    colnames(tab)[3]="Allele"
    colnames(tab)[4]="Count"
    tab["SNP"]=tab["Gene"]
    tab[tab$Allele=="All1","Allele"]="alt"
    tab[tab$Allele=="All2","Allele"]="ref"
    num<-tab %>% group_by(Gene) %>% summarise(numAlt=sum(Allele=="alt"),numRef=sum(Allele=="ref")) %>% as.data.frame()
    tab=tab[tab$Gene %in% num[num$numAlt>0 & num$numRef>0,"Gene"],]
    out=TestSNP_aod(tab,minCount=10,minSamp=0)
    return(out)
})

saveRDS(ASE_gene,"ASE_genes.RDS")