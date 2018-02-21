#!/bin/bash

set -o nounset
set -o errexit
set -o pipefail

REMOVEDUPS="/mnt/store1/TESTBISMARK/TEMP/removeN6duplicates.py"
SAMTOOLS="/home/rmgzshd/samtools-1.4.1/samtools"
BISMBAMDIR="/rdZone/live/rd00e2/BISMBAM"
POOL2DIR="/rdZone/live/rd00e2/raw_fastqs_170615_D00623_0114_BCALMRANXX_NET_Pool_2"


while read BISMBAM
        do
        PREFIX=$(echo ${BISMBAM} | sed 's/R1_001.preprocessed_bismark_bt2.bam//')
       	FQ="$PREFIX"R2_001.fastq
       	SORTEDSAM="$PREFIX"R1_001.preprocessed_bismark_bt2.sorted.sam
       	DEDUPBAM="$PREFIX"R1_001.preprocessed_bismark_bt2.sorted.dedup.bam
       	DEDUPLOG="$PREFIX"R1_001.preprocessed_bismark_bt2.sorted.dedup.log
       	
        icd $BISMBAMDIR
        iget $BISMBAM


       	icd $POOL2DIR
       	iget "$FQ".gz
       	gunzip "$FQ".gz

       	$SAMTOOLS sort -@8 -O SAM $BISMBAM |
       	grep -v "Sailfish" > $SORTEDSAM

        python $REMOVEDUPS -i $SORTEDSAM -b $FQ -o $DEDUPBAM -l $DEDUPLOG

        rm $BISMBAM
        rm $FQ
        rm $SORTEDSAM 

    done < BISMBAMS



