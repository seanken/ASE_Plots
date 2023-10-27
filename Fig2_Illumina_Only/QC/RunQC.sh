#! /bin/bash

#$ -cwd
#$ -e ErrFiles/QC.err
#$ -o ErrFiles/QC.log
#$ -q broad
#$ -P regevlab
#$ -l h_vmem=90g
#$ -l h_rt=60:00:00
#$ -l os="RedHat7"

source /broad/software/scripts/useuse

source /home/unix/ssimmons/ForPyth.sh
use Java-1.8
#use .java-jdk-1.8.0_181-x86-64
use .hdf5-1.8.16
use .picard-tools
use BEDTools

#samp=YJ89
#samp=Z93S
samp=$1
bam=/stanley/levin_asap_storage/612-eqtl/DS_ForPaper/GTExData/bams/samp_130_${samp}_130.bam
cells=/stanley/levin_asap_storage/612-eqtl/DS_ForPaper/GTExData/samp_${samp}_130/output/STARSolo/output/resultsSolo.out/GeneFull/raw/barcodes.tsv
gtf=/stanley/levin_asap/ssimmons/eQTL/Code/AlleleExpressionPipeline/version2.0/ASE_pipeline/ref/genes.gtf
prefix=/broad/hptmp/ssimmons/QTL/GTExData/QC/samp_130_${samp}_130.forQC

jarfile=/stanley/levin_dr/ssimmons/GeneralCode/CellLevel_QC/Updated/Jar/SingleCellQC.jar

echo Preprocess
#source /stanley/levin_dr/ssimmons/GeneralCode/CellLevel_QC/Updated/scripts/PrepSTARSolo.sh $bam $gtf $prefix

echo Get QC
java -jar $jarfile -c $cells -g $gtf -i ${prefix}.label.bam -o ${samp}.QC.txt -q STARSolo -v
