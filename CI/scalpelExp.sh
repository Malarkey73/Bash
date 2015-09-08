#!/bin/bash

GENOME="/mnt/store1/GATK_BUNDLE_2.8_hg19/Homo_sapiens_assembly19.fasta"
BED="/home/rmgzshd/symlinks/store1/WXS_CAPTURE_BEDS/TLX1_TLX3_SR.bed"
SCALPEL="/home/rmgzshd/scalpel-0.4.1/scalpel"

export GENOME; export BED; export SCALPEL

for BAM in *.bam
do
	PREFIX=$(echo ${BAM} | sed 's/_L001.bam//')
	$SCALPEL --export --db $PREFIX/variants.db --type all --bed $BED --covratio 0.1 --mincov 30 --ref $GENOME --format vcf > $PREFIX/$PREFIX.vcf 
done
