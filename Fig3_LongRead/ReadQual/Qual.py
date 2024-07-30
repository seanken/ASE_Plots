import pysam
def GetQual(bamfile,illumina=False):
    samfile = pysam.AlignmentFile(bamfile, "rb")
    vals=["match","insert","del","mismatch"]
    ret=[0 for i in vals]
    iter=0
    for seq in samfile:
        iter=iter+1
        if iter % 10000==0:
            print(iter)
            print(ret)
            tot=sum(ret)
            print([100*r/tot for r in ret])
            print(" ")
        cigar=seq.cigartuples
        if cigar==None:
            continue;
        ret[0]=sum([c[1] for c in cigar if c[0]==0 or c[0]==7])+ret[0]
        ret[1]=sum([c[1] for c in cigar if c[0]==1])+ret[1]
        ret[2]=sum([c[1] for c in cigar if c[0]==2])+ret[2]
        ret[3]=sum([c[1] for c in cigar if c[0]==8])+ret[3]
        if illumina:
            if not seq.has_tag("nM"):
                continue;
            nM=seq.get_tag("nM")
            ret[0]=ret[0]-nM
            ret[3]=ret[3]+nM
    return(ret)



    
