#! /bin/bash

#$ -cwd
#$ -e ErrFiles/star.$TASK_ID.err
#$ -o ErrFiles/star.$TASK_ID.log
#$ -q broad
#$ -P regevlab
#$ -l h_vmem=60g
#$ -l h_rt=60:00:00
#$ -l os="RedHat7"

#$ -t 1-2

source /broad/software/scripts/useuse

source /home/unix/ssimmons/ForPyth.sh

use BEDTools
use .bedops-2.4.14



gtf=/stanley/levin_asap/ssimmons/eQTL/Code/AlleleExpressionPipeline/version2.0/ASE_pipeline/ref/genes.gtf
SEEDFILE=samps.txt

bam=$(awk "NR==$SGE_TASK_ID" $SEEDFILE | awk '{print $1}')
workdir=$(awk "NR==$SGE_TASK_ID" $SEEDFILE | awk '{print $2}')
outfile=$(awk "NR==$SGE_TASK_ID" $SEEDFILE | awk '{print $3}')
mkdir $workdir
bed1=$workdir/UTR.3.bed
bed2=$workdir/UTR.5.bed
bed3=$workdir/Exon.bed
bed4=$workdir/Gene.bed
bed5=$workdir/Anti.bed

SEEDFILE=



echo Make Bed files
##gtf2bed in bedops
awk '{if($3=="three_prime_utr"){print $0}}' $gtf | gtf2bed > $bed1
awk '{if($3=="five_prime_utr"){print $0}}' $gtf | gtf2bed > $bed2
awk '{if($3=="exon"){print $0}}' $gtf | gtf2bed > $bed3
awk '{if($3=="gene"){print $0}}' $gtf | gtf2bed > $bed4
sed -e 's/+/-/g' -e t -e 's/-/+/g' $bed4 > $bed5

echo Tag Bam
bedtools tag -i $bam -s -f 1 -files $bed1 $bed2 $bed3 $bed4 $bed5 -labels UTR3 UTR5 Exon Gene Anti -tag GR > $workdir/ann.bam

echo Count
conda activate withPysam
python Count.Tags.py $workdir/ann.bam $outfile
