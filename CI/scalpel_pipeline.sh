#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

GENOME="/mnt/store1/GATK_BUNDLE_2.8_hg19/Homo_sapiens_assembly19.fasta"
SCALPEL="/home/rmgzshd/scalpel-0.4.1/scalpel"
TABLE_ANNOVAR="/home/rmgzshd/annovar/table_annovar.pl"
#annotations
DB="/home/rmgzshd/annovar/humandb"
#inputs


export GENOME; export SCALPEL; export TABLE_

for BAM in *.bam
do
	PREFIX=$(echo ${BAM} | sed 's/.bam//')
	BED="$PREFIX"_peaks.narrowPeak
	$SCALPEL --export --db $PREFIX/variants.db --type all --bed $BED --covratio 0.1 --mincov 30 --ref $GENOME --format vcf > $PREFIX/$PREFIX.vcf

	$TABLE_ANNOVAR $PREFIX/$PREFIX.vcf $DB \
	-buildver hg19 \
	-out $PREFIX/$PREFIX \
	-protocol refGene,1000g2015aug_all,snp138,cosmic70 \
	-operation g,f,f,f \
	-nastring . \
	-vcfinput
done