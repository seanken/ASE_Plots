library(Sierra)

tab=data.frame(Peak_file = system("ls samp_*/o*/P*/*txt",intern=T),Identifier=c("YJ89","Z93S"))
outfile="merged.peaks.txt"
#MergePeakCoordinates(tab,outfile,ncores=1)

print("Annotate")
annfile="annotated.peaks.txt"
#genome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
gtf="/stanley/levin_asap/ssimmons/eQTL/Code/AlleleExpressionPipeline/version2.0/ASE_pipeline/ref/genes.gtf"
genome<- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
AnnotatePeaksFromGTF(peak.sites.file = outfile,gtf.file = gtf,genome=genome,output.file=annfile)
