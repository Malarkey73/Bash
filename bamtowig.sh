#!/bin/bash

# data shortcuts
BAMFOLDER="/mnt/store1/rawdata/fastq/Suzana/BAM"
BEDFOLDER="/mnt/store1/rawdata/fastq/Suzana/BAM/BED"
WIGFOLDER="/mnt/store1/rawdata/fastq/Suzana/BAM/WIG"

# tool shortcuts - NB bowtie needs to be in PATH because it isn't specified in script
BEDTOOLS="/home/rmgzshd/bedtools-2.17.0/bin/bedtools"
BEDGBW="/home/rmgzshd/UCSCtools/bedGraphToBigWig"

# reference
CHROMSIZES="/home/rmgzshd/UCSCtools/hg19.chrom.sizes"

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


