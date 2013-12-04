#!/bin/bash

#shortcuts to tools
FASTQC="/home/rmgzshd/FastQC/fastqc"
# data shortcuts
FQFOLDER="/mnt/store1/rawdata/FASTQ/Sevil"

# make sure child processes see shortcuts too
export FASTQC; export SEQTK; export FQFOLDER

# 1. loop through each file in FASTQCFOLDER running FASTQC
FQS=$FQFOLDER/*.fastq.gz
for fq in $FQS
do
	$FASTQC $fq  

done
#rm $FQFOLDER/*.zip




