#! /bin/bash

#$ -cwd
#$ -e ErrFiles/star.err
#$ -o ErrFiles/star.log
#$ -q broad
#$ -P regevlab
#$ -l h_vmem=60g
#$ -l h_rt=120:00:00
#$ -l os="RedHat7"

source /broad/software/scripts/useuse

use BEDTools
use .samtools-1.15.1
use .bedops-2.4.14

bam=$1
#baits=$2
baits=/stanley/levin_asap/ssimmons/eQTL/PlotsForPaper/Tables/Test/baits.bed


dir=$(echo $bam | awk -F '/' '{print $7}')
mkdir $dir
echo $dir
cd $dir

echo Get uniquelly mapped
samtools view -H $bam > header.txt
samtools view $bam | grep -P 'NH:i:1\t'  >  temp.sam 

cat header.txt temp.sam | samtools view -Sb > temp.bam
#rm temp.sam

echo Make bam to bed
bedtools bamtobed -split -i temp.bam > reads.bed
#rm temp.bam


echo sort
sort-bed --max-mem 10G reads.bed > read.sort.bed
sort-bed --max-mem 10G $baits > baits.sort.bed
#rm reads.bed


echo Get closest
bedtools closest -a read.sort.bed -b baits.sort.bed -d -t first > close.bed

