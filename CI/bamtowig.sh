#!/bin/bash

# data shortcuts
BAMFOLDER="/mnt/store1/rawdata/fastq/Suzana/bam"
BEDFOLDER="/mnt/store1/rawdata/fastq/Suzana/bam/bed"
WIGFOLDER="/mnt/store1/rawdata/fastq/Suzana/bam/wig"


# tool shortcuts - NB bowtie needs to be in PATH because it isn't specified in script
BAMtoBED="/home/rmgzshd/bedtools-2.17.0/bin/bamToBed"
BEDGtoBW="/home/rmgzshd/UCSCtools/bedGraphToBigWig"

export BAMFOLDER; export BEDFOLDER; export WIGFOLDER; export BAMtoBED; export BEDtoBW;

 
BAMS=$BAMFOLDER/*.bam
for bam in $BAMS
do
	prefix=$(echo ${bam} | sed 's/.bam//')
	$BAMtoBED genomecov --ibam $bam -bg > $prefix.bedGraph  

done
rm $FQFOLDER/*.zip


