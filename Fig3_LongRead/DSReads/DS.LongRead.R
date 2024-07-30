library(ggplot2);library(dplyr);library(tidyr);library(purrr)

print("Old")
lst=system("ls /stanley/levin_asap_storage/612-eqtl/SingleNuc_data/UpdatedAIData/LongRead/Jan5Analysis/LongRead/DS/s*/o*/A*/*txt",intern=T)
tab=read.table("/stanley/levin_asap_storage/612-eqtl/SingleNuc_data/UpdatedAIData/LongRead/Jan5Analysis/LongRead/DS/bams.txt")
tab["Num"]=1:10
tab["Samp"]=sub("GTEX-","",tab[,2])
tab["Type"]=c("MAS-Seq","MAS-Seq","GridION","GridION","PacBio","PacBio",rep("PromethION",4))

tab[8,"Type"]="ProLong"
tab[10,"Type"]="ProLong"

tab=tab[c(3:7,9),]

dat=data.frame(fil=lst)

dat["DS"]=map_chr(dat[,"fil"],function(x) strsplit(x,"_")[[1]][5])

dat["Num"]=map_chr(dat[,"fil"],function(x) strsplit(strsplit(x,"_")[[1]][6],"/",fixed=T)[[1]][1])
dat["DS"]=as.numeric(dat[,"DS"]);dat["Num"]=as.numeric(dat[,"Num"])
#dat[dat$DS!=1,]
#dat[dat$Num!=8 & dat$Num!=10,]

comb=inner_join(dat,tab)
comb=comb[,c("fil","DS","Samp","Type")]


print("Revio")
lst=system("ls /stanley/levin_asap_storage/ssimmons/Revio/MA*/DS/s*/o*/A*/*txt",intern=T)
dat=data.frame("fil"=lst)
dat["DS"]=map_chr(dat[,"fil"],function(x) strsplit(x,"_")[[1]][4])
dat["Samp"]=rep(c("Z93S","YJ89"),9)
dat["Type"]="MAS-Seq+Revio"
comb=rbind(comb,dat)
lst=system("ls /stanley/levin_asap_storage/ssimmons/Revio/MA*/Re*/s*/o*/A*/*txt",intern=T)
dat=data.frame("fil"=lst)
dat["DS"]=1
dat["Samp"]=c("Z93S","YJ89")
dat["Type"]="MAS-Seq+Revio"
comb=rbind(comb,dat)

print("MAS-seq")
lst=system("ls /broad/hptmp/ssimmons/QTL/MASSeq/Skera/Pipeline/DS/s*/o*/A*/*txt",intern=T)
dat=data.frame("fil"=lst)
dat["DS"]=map_chr(dat[,"fil"],function(x) strsplit(x,"_")[[1]][4])
dat["Samp"]=rep(c("Z93S","YJ89"),8)
dat["Type"]="MAS-Seq"
comb=rbind(comb,dat)
lst=system("ls /stanley/levin_asap_storage/ssimmons/MA*/P*/samp_*/o*/A*/*txt",intern=T)
dat=data.frame("fil"=lst)
dat["DS"]=1
dat["Samp"]=c("Z93S","YJ89")
dat["Type"]="MAS-Seq"
comb=rbind(comb,dat)

print("Illumina")

dirs=c("/stanley/levin_asap_storage/612-eqtl/DS_ForPaper/GTExData/samp_YJ89_130/output","/stanley/levin_asap_storage/612-eqtl/DS_ForPaper/GTExData/samp_Z93S_130/output")
names(dirs)=c("YJ89","Z93S")

fils=map(dirs,function(x){
    fil=map_chr(1:9,function(y) paste(x,"/Downsample/counts.",y,".txt",sep=""))
    fil=c(fil,paste(x,"/AlleleCounts/counts.txt",sep=""))
    DS=.1*1:10
    dat=data.frame(fil,DS)
})

for(i in names(fils)){fils[[i]]["Samp"]=i}

fils=do.call(rbind,fils)
fils["Type"]="Illumina"

#fils["Num"]=c(rep(11,10),rep(12,10))

fils=rbind(comb[,colnames(fils)],fils)

fils["Number"]=map_dbl(fils[,"fil"],function(x){print(x);tab=read.table(x);tab=tab[tab[,3]!="Ambig",];sum(tab[,4])})

#fils=fils[fils$Num!=8 & fils$Num!=10,]
fils["DS"]=as.numeric(fils[,"DS"])


fils["NumReads"]=0
fils["summary"]=sub("AlleleCounts/counts.txt","FC_exon/gene_assigned.summary",fils[,"fil"],fixed=T)
fils[fils$Type!="Illumina","NumReads"]=map_dbl(fils[fils$Type!="Illumina","summary"],function(x){tab=read.table(x,header=T);print(head(tab));sum(tab[,2])})

fils["SumSTARSolo"]=map_chr(fils[,"fil"],function(x) paste0(c(strsplit(x,"/",fixed=T)[[1]][1:8],"STARSolo/output/resultsSolo.out/GeneFull/Summary.csv"),collapse="/") )

fils[fils$Type=="Illumina","NumReads"]=map_dbl(fils[fils$Type=="Illumina","SumSTARSolo"],function(x) read.table(x,sep=",")[1,2])

fils[fils$Type=="Illumina","NumReads"]=apply(fils[fils$Type=="Illumina",c("NumReads","DS")],1,function(x) x["NumReads"]*x["DS"])

saveRDS(fils,"fils.withlong.RDS")
#fils["Type"]=factor(fils[,"Type"],unique(fils[,"Type"])[c(5,2,4,3,1)])

p=ggplot(fils,aes(x=NumReads,y=Number,color=Type))+geom_point()+geom_line()+facet_wrap(~Samp)+scale_x_log10()+scale_y_log10()+annotation_logticks(sides = "lb")+xlab("Number of reads")+ylab("Number of phased UMIs")
ggsave("DSLongReads.pdf",p,width=14)
