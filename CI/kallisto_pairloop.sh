#!/bin/bash
KALLISTO="/home/rmgzshd/kallisto/kallisto"
TRANSCRIPTS="/mnt/store1/Pathania/Hs.GRCh38.idx"


export KALLISTO; export TRANSCRIPTS;

for fq in *_1.fq.gz
do
	PREFIX=$(echo ${fq} | sed 's/_1.fastq.gz//')
	$KALLISTO quant -i $TRANSCRIPTS -o $PREFIX	-b 100 "$PREFIX"_1.fastq.gz "$PREFIX"_2.fastq.gz
done