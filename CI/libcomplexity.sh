#!/bin/bash

# data shortcuts
BAMFOLDER="/mnt/store1/rawdata/fastq/SA1SA2/BAM"

#tools shortcuts
C_CURVE="/home/rmgzshd/preseq-0.0.4/c_curve"
LC_EXTRAP="/home/rmgzshd/preseq-0.0.4/lc_extrap"
PPCQTOOLS="/home/rmgzshd/phantompeakqualtools/run_spp.R"

export BAMFOLDER; export C_CURVE; export LC_EXTRAP

mkdir libcomplexity 

for bam in $BAMFOLDER/*.bam
do
	prefix=$(echo ${bam} | sed 's/.bam//')
	$C_CURVE -B $bam -o $prefix.ccurve.txt
	$LC_EXTRAP -B $bam -o $prefix.lcextrap.txt
	Rscript $PPCQTOOLS -c=$bam -savp -out=$prefix.ppcq.txt
# Rscript /home/rmgzshd/phantompeakqualtools/run_spp.R -c=SRR445739.sorted.bam -savp -out=SRR445739.sorted.ppcq.txt
done
mkdir -p libcomplexity
mv $BAMFOLDER/*.ccurve.txt
mv $BAMFOLDER/*.lcextrap.txt
mv $BAMFOLDER/*.ppcq.txt
