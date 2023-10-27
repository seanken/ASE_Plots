library(scAlleleExpression)
library(Seurat)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(purrr)



fil_phase="/stanley/levin_asap_storage/612-eqtl/DS_ForPaper/GTExData/samp_Z93S_130/output/AlleleCounts/counts.txt"
fil_vcf="/stanley/levin_asap_storage/612-eqtl/DS_ForPaper/GTExData/samp_Z93S_130/output/SNPLevelAlleleCounts/counts.txt"
fil_novcf=""

fils=c(fil_phase,fil_vcf,fil_novcf)

snpToGene=read.table("/stanley/levin_asap_storage/612-eqtl/DS_ForPaper/GTExData/samp_Z93S_130/output/SNPLevelCounts/comb.bed")
dat_phase=read.table(fil_phase)
dat_phase=dat_phase[dat_phase$V3!="Ambig",]

dat_phase <- dat_phase %>% group_by(V2) %>% summarise(Num=sum(V4))

dat_vcf=read.table(fil_vcf)

snpToGene=snpToGene[,1:2]
colnames(snpToGene)=c("SNP","Gene")

colnames(dat_vcf)[2]="SNPName"

snpToGene["SNPName"]=map_chr(snpToGene[,"SNP"],function(x){x=sub("chr","",x);s=strsplit(x,"_")[[1]];paste(s[1],(as.numeric(s[2])+1),sep=":")})

snpToGene=snpToGene[!duplicated(snpToGene[,"SNPName"]),]

dat_vcf=dat_vcf[dat_vcf$V3!="Ambig",]

dat_vcf <- dat_vcf %>% group_by(SNPName) %>% summarise(Num=sum(V4)) %>% as.data.frame()

dat_vcf=inner_join(dat_vcf,snpToGene)