#!/bin/bash


#$ -cwd
#$ -q broad
#$ -P regevlab
#$ -l h_vmem=90g
#$ -e ErrFile/err.$TASK_ID.err
#$ -o ErrFile/log.$TASK_ID.out
#$ -l os=RedHat7
#$ -l h_rt=400:00:00

#$ -t 1-8

source /broad/software/scripts/useuse

use UGER
source ~/ForPyth.sh

#conda activate Milo 
conda activate AlleleDownstream2

Rscript CompareDAESC.only.R $SGE_TASK_ID
#Rscript CompareDAESC.R
