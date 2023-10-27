#!/bin/bash


#$ -cwd
#$ -q broad
#$ -P regevlab
#$ -l h_vmem=90g
#$ -e ErrFile/err.err
#$ -o ErrFile/log.out
#$ -l os=RedHat7
#$ -l h_rt=140:00:00

source /broad/software/scripts/useuse

use UGER
source ~/ForPyth.sh

#conda activate Milo 
conda activate AlleleDownstream2

Rscript TestDAESC.R
