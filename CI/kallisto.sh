#!/bin/bash

INDEX="/mnt/store1/KALLISTO_INDEX/Mus_musculus.GRCm38.ERCC.eGFP.rel79.cdna.all.fa.idx"

R1="_R1_001.fastq.gz"
R2="_R2_001.fastq.gz"

for fq in *R1_001.fastq.gz
do
	PREFIX=$(echo ${fq} | sed 's/_R1_001.fastq.gz//')
	kallisto quant -i $INDEX -o $PREFIX -t 20 $PREFIX$R1 $PREFIX$R2
done

