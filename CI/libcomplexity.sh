#!/bin/bash

# data shortcuts
BAMFOLDER="/mnt/store1/rawdata/FASTQ/SA1SA2/BAM"

#tools shortcuts
C_CURVE="/home/rmgzshd/preseq-0.0.4/c_curve"
LC_EXTRAP="/home/rmgzshd/preseq-0.0.4/lc_extrap"
PPCQTOOLS="/home/rmgzshd/phantompeakqualtools/run_spp.R"
CHIPQCPLOTS="/home/rmgzshd/RScripts/ChIPQCplots.R"

export BAMFOLDER; export C_CURVE; export LC_EXTRAP; export CHIPQCPLOTS


# we compute the cross correlation, the lib complexity and estimate
# extra sequencing required and ... 

for bam in $BAMFOLDER/*.sorted.bam
do
	prefix=$(echo ${bam} | sed 's/.sorted.bam//')
	$C_CURVE -B $bam -o $prefix.ccurve.txt
	$LC_EXTRAP -B $bam -o $prefix.lcextrap.txt
	Rscript $PPCQTOOLS -c=$bam -savp=$prefix.crosscor.pdf -out=$prefix.ppcq.txt
done

# the RScript itself should recognise the right files 
# and loop through converting them to plots (req ggplot2)
Rscript $CHIPQCPLOTS

#rename 's/.sorted.pdf$/.crosscor.pdf/' *.sorted.pdf
mkdir -p libcomplexity
mv $BAMFOLDER/*.pdf libcomplexity
rm $BAMFOLDER/*.ccurve.txt
rm $BAMFOLDER/*.lcextrap.txt
rm $BAMFOLDER/*.ppcq.txt
