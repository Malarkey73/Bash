#!/bin/bash

GENOME="/mnt/store1/GATK_BUNDLE_2.8_hg19/Homo_sapiens_assembly19.fasta"
SAMTOOLS="/home/rmgzshd/samtools-1.2/samtools"
BWA="/home/rmgzshd/bwa/bwa"

export GENOME; export SAMTOOLS; export BWA

# in order for this script to work the samples must bedescribed in a specific way:
# e.g. NB1.N.L2.R1.fastq.gz , sample, type, lane, pair then suffix descriptor fastq.gz
# the first 3 go into the @RG descriptor of the header that is used within GATK software
for fq in *R1.fastq.gz
do
	echo $fq
    # get the sample name prefix
   	PREFIX=$(echo ${fq} | sed 's/.R1.fastq.gz//')
	foo=${PREFIX##*.}
	$BWA mem -M -t 20 -R "@RG\tID:$PREFIX\tSM:$foo\tPL:ILLUMINA" $GENOME $PREFIX.R1.fastq.gz $PREFIX.R2.fastq.gz |
	$SAMTOOLS view -Shu - | 
	$SAMTOOLS sort -m2G -@8 - $PREFIX
done


