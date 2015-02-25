#!/bin/bash

# summary files and paths
SUMMARYFILES=`find . -print | grep -i 'summary.txt'`
FIRST=true

for SUMMARY in $SUMMARYFILES
do
  	#if this is the first sample then put a header row in first
        if [[ $FIRST = true ]]; then
                printf "SAMPLE\t" > fastqc.results
                cut -f 2 $SUMMARY | tr "\n" "\t" >> fastqc.results
                printf "\n"  >>fastqc.results
                FIRST=false
        fi

	cut -f 3 $SUMMARY | head -1 | tr "\n" "\t" >> fastqc.results
        cut -f 1 $SUMMARY | tr "\n" "\t" >> fastqc.results
        printf "\n" >> fastqc.results
done