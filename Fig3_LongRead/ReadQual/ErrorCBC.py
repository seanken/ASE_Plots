import pysam
import sys

bamfile=sys.argv[1]

samfile = pysam.AlignmentFile(bamfile, "rb")

val=0
tot=0
for seq in samfile:
    cbc_uncor=seq.get_tag("CR")
    if not seq.has_tag("CB"):
        #val=val+1;
        continue
    tot=tot+1
    cbc_cor=seq.get_tag("CB")[0:16]
    if cbc_uncor==cbc_cor:
        val=val+1
    if tot>10000000:
        break;
print(val)
print(tot)
rat=val/tot
print(rat)