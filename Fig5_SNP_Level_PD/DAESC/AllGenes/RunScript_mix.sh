#!/bin/bash


#$ -cwd
#$ -q broad
#$ -P regevlab
#$ -l h_vmem=90g
#$ -e ErrFile/err.mix.$TASK_ID.err
#$ -o ErrFile/log.mix.$TASK_ID.out
#$ -l os=RedHat7
#$ -l h_rt=140:00:00

#$ -t 1-1000
#$ -tc 100

source /broad/software/scripts/useuse

use UGER
source ~/ForPyth.sh

#conda activate Milo 
conda activate AlleleDownstream2

Rscript DAESC.Gene.level.R $SGE_TASK_ID
