import pysam
import sys
import pandas as pd
import random
import numpy as np

args=sys.argv
bamfile=args[1]
outfile=args[2]
samfile=pysam.AlignmentFile(bamfile, "rb") 
savfil=open(outfile,"w")
iter=0
for read in samfile:
    iter=iter+1
    if iter % 100000==0:
        print(iter)
    ran=random.uniform(a=0,b=1)
    if ran>.01:
        continue;
    if read.has_tag("SA") or not read.has_tag("tp"):
        continue;
    #suppmaps=read.get_tag("SA")
    isprim=read.get_tag("tp")
    if isprim!="P":
        continue;
    cigars=read.cigartuples;
    vals=[c[1] for c in cigars if c[0]==0 or c[0]==7]
    readlen=np.sum(vals)
    savfil.write(str(readlen)+"\n")
samfile.close()
savfil.close()
    



    

