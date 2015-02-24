#!/bin/bash

touch fastq.results
echo "BasicStatistics/tPerBaseSequenceQuality/tPerSequenceQualityScores"

SUMMARYFILES=`find . -print | grep -i 'summary.txt'`
for SUMMARY in $SUMMARYFILES
do
	cat $SUMMARY | 
done