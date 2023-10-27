use .bedops-2.4.14
use BEDTools

vcf=/stanley/levin_asap_storage/612-eqtl/eQTL_testing/GTEx/geno.no.chr.vcf.gz
fasta=/stanley/levin_dr/ssimmons/TRN/Intron/Intron/fasta/genome.fa

echo Extract bed
zcat $vcf > vcf2bed > geno.bed

echo Extract seq
bedtools getfasta -fi $fasta -bed geno.bed -bedOut > out.bed

echo Done
