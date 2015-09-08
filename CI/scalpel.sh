#!/bin/bash

GENOME="/mnt/store1/GATK_BUNDLE_2.8_hg19/Homo_sapiens_assembly19.fasta"
BED="/home/rmgzshd/symlinks/store1/WXS_CAPTURE_BEDS/TLX1_TLX3_SR.bed"
SCALPEL="/home/rmgzshd/scalpel-0.4.1/scalpel"

export GENOME; export BED; export SCALPEL

for BAM in *.bam
do
	PREFIX=$(echo ${BAM} | sed 's/_L001.bam//')
	$SCALPEL --single --bam $BAM --bed $BED --numprocs 20 --ref $GENOME --dir ./"$PREFIX"
done
