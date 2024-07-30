#! /bin/bash

#$ -cwd
#$ -e ErrFiles/err.err
#$ -o ErrFiles/out.log
#$ -q broad
#$ -P regevlab
#$ -l h_vmem=60g
#$ -l h_rt=120:00:00
#$ -l os="RedHat7"


#$ -t 1-10
source /broad/software/scripts/useuse

source ~/ForPyth.sh

SEEDFILE=samps.txt



use BEDTools
use .samtools-1.8
#use .java-jdk-1.8.0_181-x86-64
use Java-1.8

conda activate Milo

SEEDFILE=samples.txt

samp=$(awk "NR==$SGE_TASK_ID" $SEEDFILE | awk '{print $2')
bam=$(awk "NR==$SGE_TASK_ID" $SEEDFILE | awk '{print $1}')
