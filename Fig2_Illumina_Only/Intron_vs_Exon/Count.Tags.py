import pysam
import sys
import pandas as pd

args=sys.argv
bamfile=args[1]
outfile=args[2]
samfile=pysam.AlignmentFile(bamfile, "rb")
totCounts={}
snpCounts={} 
for read in samfile:
    if not read.has_tag("NH"):
        continue;
    num=read.get_tag("NH")
    if num!=1:
        continue;
    tagVal="Intergenic"
    iterat=0
    if read.has_tag("GR"):
        iterat=iterat+1
        if iterat % 100000==0:
            print(iterat)
        gr=read.get_tag("GR").split(";")
        
        if "Anti" in gr:
            tagVal="Antisense"
        if "Gene" in gr:
            tagVal="Intronic"
        if "Exon" in gr:
            tagVal="Exonic"
        if "UTR3" in gr:
            tagVal="3' UTR"
        if "UTR5" in gr:
            tagVal="5' UTR"
        if "UTR3" in gr and "UTR5" in gr:
            tagVal="UTR_both"
        
    
    totCounts[tagVal]=totCounts.get(tagVal,0)+1
    if read.has_tag("vW"):
        snpCounts[tagVal]=snpCounts.get(tagVal,0)+1
samfile.close()
#totRead=[totCounts.get(tagVal,0) for tagVal in totCounts.keys()]
#snpRead=[snpCounts.get(tagVal,0) for tagVal in totCounts.keys()]
#region=[tagVal for tagVal in totCounts.keys()]

dat=[[totCounts.get(tagVal,0),snpCounts.get(tagVal,0),tagVal] for tagVal in totCounts.keys()]
dat=pd.DataFrame(dat,columns=["totCounts","snpCounts","region"])
print(dat)
dat.to_csv(outfile)
