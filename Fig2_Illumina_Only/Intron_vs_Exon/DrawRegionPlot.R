library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

fil1="out.YJ89.txt"
fil2="out.Z93S.txt"
out1=read.csv(fil1)
out2=read.csv(fil2)
out1["Sample"]="YJ89"
out2["Sample"]="Z93S"

dat=rbind(out1,out2)

dat=dat[dat$region!="UTR_both",]

dat["PercentWithSNP"]=100*dat[,"snpCounts"]/dat[,"totCounts"]
p=ggplot(dat,aes(x=region,y=PercentWithSNP,fill=Sample))+geom_bar(position="dodge",stat="identity")+coord_flip()+xlab("Region")+ylab("% Uniquelly Mapped Reads Overlapping a Het SNP")
ggsave("RegionBarPlot.pdf",p,width=10)
