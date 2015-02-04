#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

#tools
SAMTOOLS="/home/rmgzshd/samtools/samtools"
VARSCAN="/home/rmgzshd/VarScan/VarScan.v2.3.7.jar"

#places
GENOME="/mnt/store1/GATK_BUNDLE_2.8_hg19/HPV_BUNDLE/GRCh37-lite-+-HPV_Redux-build.fa"
CAPTURE_REGION="/mnt/store1/WXS_CAPTURE_BEDS/SeqCap_EZ_Exome_v3_primary.bed"
DATAFOLDER="/mnt/store1/CESC_WXS"

#args
TUMORBAM=`echo $1`
REFBAM=`echo $2`
PREFIX=`echo $TUMORBAM | sed 's/.bam//' `

export GENOME; export CAPTURE_REGION; export DATAFOLDER; export TUMORBAM; export REFBAM; export PREFIX

hostname
date
mkdir -p $DATAFOLDER/pileups

MPILEUPDATA="$SAMTOOLS mpileup -A -B -C 50 -d 10000 -m 3 -F 0.0002 -q 20 -Q 20 \
-f $GENOME \
-l $CAPTURE_REGION \
$DATAFOLDER/$TUMORBAM \
$DATAFOLDER/$REFBAM > \
$DATAFOLDER/pileups/$PREFIX.mpileup"

java -d64 -Xmx4g -jar $VARSCAN somatic <\($MPILEUPDATA\) $PREFIX \
--mpileup 1 --min-coverage 8 --min-coverage-normal 10 --min-coverage-tumor 6 \
--min-var-freq 0.01 --min-freq-for-hom 0.75 --normal-purity 1 --p-value 0.99 \
--somatic-p-value 0.05 --tumor-purity 0.5 --strand-filter 0 
