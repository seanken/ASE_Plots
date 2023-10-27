#!/bin/bash


#$ -cwd
#$ -q broad
#$ -P regevlab
#$ -l h_vmem=40g
#$ -e ErrFile/err.$TASK_ID.err
#$ -o ErrFile/log.$TASK_ID.out
#$ -l os=RedHat7
#$ -l h_rt=40:00:00

#$ -t 1-609
#$ -tc 30

source /broad/software/scripts/useuse

use UGER
source ~/ForPyth.sh

conda activate AlleleDownstream2 

Rscript Trans.eQTL.R $SGE_TASK_ID
