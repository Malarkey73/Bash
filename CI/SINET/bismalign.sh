#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

BISMARK="/home/rmgzshd/Bismark/bismark"
GENOMEFOLDER="/mnt/store1/hg38"

for FQ in *.preprocessed.fq
	do
	PREFIX=$(echo ${FQ} | sed 's/_R1_001.preprocessed.fq//')
	$BISMARK --bowtie2 --sam --parallel 2 --basename $PREFIX $GENOMEFOLDER $FQ


	done