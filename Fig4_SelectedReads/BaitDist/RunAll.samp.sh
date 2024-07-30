for bam in  /stanley/levin_asap_storage/6*/n*21*/CellRanger/*SM*/s1/outs/*bam; do qsub CheckDist.samp.sh $bam; done
