library(scAlleleExpression)
library(Seurat)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(purrr)



fil_phase="/stanley/levin_asap_storage/612-eqtl/DS_ForPaper/GTExData/samp_Z93S_130/output/AlleleCounts/counts.txt"
fil_vcf="/stanley/levin_asap_storage/612-eqtl/DS_ForPaper/GTExData/samp_Z93S_130/output/SNPLevelAlleleCounts/counts.txt"
fil_novcf="/stanley/levin_asap_storage/612-eqtl/SNPsCalling/CellSNP/RunPipe/samp_Z93S/output/SNPLevelAlleleCounts/counts.txt"

fils=c(fil_phase,fil_vcf,fil_novcf)

snpToGene=read.table("/stanley/levin_asap_storage/612-eqtl/DS_ForPaper/GTExData/samp_Z93S_130/output/SNPLevelCounts/comb.bed")
dat_phase=read.table(fil_phase)
dat_phase=dat_phase[dat_phase$V3!="Ambig",]

dat_phase <- dat_phase %>% group_by(V2) %>% summarise(Num=sum(V4)) %>% data.frame()

dat_phase["Type"]="Phase"

dat_vcf=read.table(fil_vcf)

snpToGene=snpToGene[,1:2]
colnames(snpToGene)=c("SNP","Gene")

colnames(dat_vcf)[2]="SNPName"

snpToGene["SNPName"]=map_chr(snpToGene[,"SNP"],function(x){x=sub("chr","",x);s=strsplit(x,"_")[[1]];paste(s[1],(as.numeric(s[2])-1),sep=":")})

#snpToGene=snpToGene[!duplicated(snpToGene[,"SNPName"]),]

dat_vcf=dat_vcf[dat_vcf$V3!="Ambig",]

dat_vcf <- dat_vcf %>% group_by(SNPName) %>% summarise(Num=sum(V4)) %>% as.data.frame()

dat_vcf=inner_join(dat_vcf,snpToGene)

dat_vcf["Type"]="NoPhase"

print("No VCF!")
dat_novcf=read.table(fil_novcf)

#snpToGene=snpToGene[,1:2]
#colnames(snpToGene)=c("SNP","Gene")

colnames(dat_novcf)[2]="SNPName"

#snpToGene["SNPName"]=map_chr(snpToGene[,"SNP"],function(x){x=sub("chr","",x);s=strsplit(x,"_")[[1]];paste(s[1],(as.numeric(s[2])-1),sep=":")})

#snpToGene=snpToGene[!duplicated(snpToGene[,"SNPName"]),]

dat_novcf=dat_novcf[dat_novcf$V3!="Ambig",]

dat_novcf <- dat_novcf %>% group_by(SNPName) %>% summarise(Num=sum(V4)) %>% as.data.frame()

dat_novcf=inner_join(dat_novcf,snpToGene)


dat_novcf["Type"]="NoVCF"

ret=list()
ret[["Phase"]]=dat_phase
ret[["NoPhase"]]=dat_vcf
ret[["NoVCF"]]=dat_novcf

saveRDS(ret,"ret.RDS")

ret=ret %>% map(function(x) x %>% group_by(Gene,Type) %>% summarise(Num=max(Num)) %>% as.data.frame())


dat=do.call(rbind,ret) %>% spread(Type,Num,fill=0) %>% gather(Type,Count,-Gene)

dat["Title"]="VCF, phased"
dat[dat$Type=="NoPhase","Title"]="VCF, not phased"
dat[dat$Type=="NoVCF","Title"]="no VCF"
dat["Title"]=factor(dat[,"Title"],unique(dat[,"Title"])[c(3,1,2)])
ggplot(dat,aes(y=Count+1,x=Title,fill=Title))+theme(legend.position="none")+scale_y_log10()+geom_violin(scale="width")+ylab("Number of phased UMIs plus 1")+xlab("")+annotation_logticks(sides = "l")
ggsave("VCFType.pdf")