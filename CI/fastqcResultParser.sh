#!/bin/bash

# given a folder full of fastqc results from multiple samples this script will
# make a simpel summary table of all those results in a single file
# you might use it a s part of a QC process with the similar flagstatResultParser.sh
# which is for BAM file stats.

# summary files and paths
SUMMARYFILES=`find . -print | grep -i 'summary.txt'`
FIRST=true

for SUMMARY in $SUMMARYFILES
do
  		#if this is the first sample then print a header row
        if [[ $FIRST = true ]]; then
                printf "SAMPLE\t" > fastqc.summary
                cut -f 2 $SUMMARY | tr "\n" "\t" >> fastqc.summary
                printf "\n"  >>fastqc.summary
                FIRST=false
        fi
    # the name and summary resulst for each sample 1 row each
	cut -f 3 $SUMMARY | head -1 | tr "\n" "\t" >> fastqc.summary
    cut -f 1 $SUMMARY | tr "\n" "\t" >> fastqc.summary
    printf "\n" >> fastqc.summary
done