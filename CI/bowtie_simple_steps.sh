#!/bin/bash

# data shortcuts
FQFOLDER="/mnt/store1/rawdata/fastq/Suzana"
BAMFOLDER="/mnt/store1/rawdata/fastq/Suzana/bam"
# reference shortcuts
GENOME="/mnt/store1/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome"
GENES="/mnt/store1/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf"
# tool shortcuts
SEQTK="/home/rmgzshd/seqtk/seqtk"
BOWTIE2="/home/rmgzshd/bowtie2-2.1.0/bowtie2"
SAMTOOLS="/home/rmgzshd/samtools/samtools"

export FQFOLDER; export BAMFOLDER; export GENOME; export GENES; export SEQTK; export BOWTIE2; export SAMTOOLS

# trim
$SEQTK trimfq $FQFOLDER/SC1_A6_R1_001.fq.gz > $FQFOLDER/SC1_A6_R1_001.trim.fq
$SEQTK trimfq $FQFOLDER/SC1_A6_R2_001.fq.gz > $FQFOLDER/SC1_A6_R2_001.trim.fq
# align
mkdir -p $BAMFOLDER
$BOWTIE2 -p 24 -x $GENOME -1 $FQFOLDER/SC1_A6_R1_001.trim.fq -2 $FQFOLDER/SC1_A6_R2_001.trim.fq -S $BAMFOLDER/SC1_A6.sam
# sam to bam
$SAMTOOLS view -bS $BAMFOLDER/SC1_A6.sam -o $BAMFOLDER/SC1_A6.bam	
# sort
$SAMTOOLS sort $BAMFOLDER/SC1_A6.bam $BAMFOLDER/SC1_A6.sorted
# tidy
rm $BAMFOLDER/*.sam
rm $BAMFOLDER/SC1_A6.bam


