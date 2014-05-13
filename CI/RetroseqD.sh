#!/bin/bash

# data shortcuts
TARGETBAM="/mnt/store/cghub/HNSC/62d09912-e98d-4611-be72-22cb4cd8b6b6/TCGA-CN-5367-01A-01D-1431_120420_SN1120_0134_AC0J8YACXX_s_2_rg.sorted.bam"

RETROSEQ="/home/rmgzshd/RetroSeq/bin/retroseq.pl"

export TARGETBAM; export RETROSEQ;

prefix=$(echo ${TARGETBAM} | sed 's/.sorted.bam//')

STARTTIME=$(date +%s)

$RETROSEQ -discover -bam $TARGETBAM -output $prefix.candidates.tab -refTEs ref_types.tab -eref probes.tab -align

ENDTIME=$(date +%s)
dt=$(($ENDTIME - $STARTTIME))
ds=$((dt % 60))
dm=$(((dt / 60) % 60))
dh=$((dt / 3600))
printf '\n\t %d:%02d:%02d \n' $dh $dm $ds

