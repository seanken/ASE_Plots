#! /bin/bash

#$ -cwd
#$ -e ErrFiles/err.err
#$ -o ErrFiles/out.log
#$ -q broad
#$ -P regevlab
#$ -l h_vmem=60g
#$ -l h_rt=120:00:00
#$ -l os="RedHat7"

source /broad/software/scripts/useuse

source ~/ForPyth.sh

SEEDFILE=samps.txt



use BEDTools
use .samtools-1.8
#use .java-jdk-1.8.0_181-x86-64
use Java-1.8

export PATH=$PATH:/stanley/sheng-lab/ssimmons/Grin2/SingleNuc/IsoformUsage/Sierra/New/regtools/regtools/build
conda activate Milo

nextflow=/stanley/levin_asap/ssimmons/eQTL/Code/AlleleExpressionPipeline/version2.0/ASE_pipeline/nextflow
#pipeline=/stanley/sheng-lab/ssimmons/Grin2/SingleNuc/IsoformUsage/Sierra/New/Code/GetPeaks.nf
pipeline=/stanley/levin_asap/ssimmons/eQTL/Code/AlleleExpressionPipeline/RunSierra/GetPeaks.nf

#SEEDFILE=samps.txt
samp=$1
indir=/broad/hptmp/ssimmons/QTL/GTExData/samp_130_${samp}_130/output
#samp=$(awk "NR==$SGE_TASK_ID" $SEEDFILE | awk '{print $1}')
#indir=$(awk "NR==$SGE_TASK_ID" $SEEDFILE | awk '{print $3}')
bam=$indir/STARSolo/output/resultsAligned.sortedByCoord.out.bam
cells=$indir/STARSolo/output/resultsSolo.out/GeneFull/filtered/barcodes.tsv
gtf=/stanley/levin_asap/ssimmons/eQTL/Code/AlleleExpressionPipeline/version2.0/ASE_pipeline/ref/genes.gtf

dir=samp_$samp
mkdir $dir
cd $dir

outdir=output

$nextflow $pipeline --bam $bam --gtf $gtf --cells $cells --outdir $outdir

