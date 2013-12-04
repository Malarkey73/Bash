#!/bin/bash

# data shortcuts
BAMFOLDER="/mnt/store1/rawdata/FASTQ/SA1SA2/BAM"
BEDFOLDER="/mnt/store1/rawdata/FASTQ/SA1SA2/BAM/BED"
WIGFOLDER="/mnt/store1/rawdata/FASTQ/SA1SA2/BAM/WIG"

# tool shortcuts
BEDTOOLS="/home/rmgzshd/bedtools-2.17.0/bin/bedtools"
BEDGBW="/home/rmgzshd/UCSCtools/bedGraphToBigWig"

# reference
CHROMSIZES="/home/rmgzshd/UCSCtools/mm10.chrom.sizes"

export BAMFOLDER; export BEDFOLDER; export WIGFOLDER; export BEDTOOLS; export BEDGBW;

 
BAMS=$BAMFOLDER/*.bam
for bam in $BAMS
do
	prefix=$(echo ${bam} | sed 's/.bam//')

	$BEDTOOLS genomecov -ibam $bam -bg > $prefix.bedGraph
	$BEDGBW $prefix.bedGraph $CHROMSIZES $prefix.bigWig

done
mkdir -p $BEDFOLDER
mkdir -p $WIGFOLDER
mv *bedGraph $BEDFOLDER
mv *bigWig $WIGFOLDER


