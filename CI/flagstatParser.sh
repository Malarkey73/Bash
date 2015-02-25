#!/bin/bash

# given a folder full of flagstat results saved in .flagtat files
# this parses them all into a single table with each row of results per sample

FLAGFILES=`ls *.flagstat`
FIRST=true
HEADER="sample\ttotal\tsecondary\tsupplementary\tduplicates\tmapped\tpaired\tread1\tread2\tproperly_paired\twith_itself\tsingletons\tmapped_diff_chr\tmapped_diff_chr_mapQ>5\n"

for FLAGFILE in $FLAGFILES
do
  		#if this is the first sample then print a header row
        if [[ $FIRST = true ]]; then
                printf $HEADER > flagstat.results
        fi
    # the name and summary resulst for each sample 1 row each
	printf $FLAGFILE >> flagstat.results
    cut -f 1 -d " " $FLAGFILE | tr "\n" "\t" >> flagstat.results
    printf "\n" >> flagstat.results
done
