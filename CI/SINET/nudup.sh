#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

icd /rdZone/live/rd00e2/BISMFQ

PYTHON27="/usr/local/bin/python2.7"
NUDUP="/mnt/store1/TESTBISMARK/nudup.py"


while read COL1
do

	iget $COL1
	PREFIX=$(echo ${COL1} | sed 's/.fq//')
	BAM="$PREFIX"_bismark_bt.bam

	$PYTHON27 $NUDUP -f $FQ -o $PREFIX $BAM

done < GETFQFROMIRODS

