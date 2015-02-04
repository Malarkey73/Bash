#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

#shortcuts to tools
SAMTOOLS="/home/rmgzshd/samtools/samtools"

# shortcuts to annotation
GENOME="/mnt/store1/GATK_BUNDLE_2.8_hg19/HPV_BUNDLE/GRCh37-lite-+-HPV_Redux-build.fa"
CAPTURE_REGION="/mnt/store1/WXS_CAPTURE_BEDS/SeqCap_EZ_Exome_v3_primary.bed"
DATAFOLDER="/mnt/store1/CESC_WXS"
TUMORBAM=`echo $1`
REFBAM=`echo $2`
PREFIX=`echo $TUMORBAM | sed 's/.bam//' `
export GENOME; export CAPTURE_REGION; export DATAFOLDER; export TUMORBAM; export REFBAM; export PREFIX




hostname
date
mkdir -p $DATAFOLDER/pileups

$SAMTOOLS mpileup -A -B -C 50 -d 10000 -m 3 -F 0.0002 -q 20 -Q 20 \
-f $GENOME \
-l $CAPTURE_REGION \
$DATAFOLDER/$TUMORBAM \
$DATAFOLDER/$REFBAM > \
$DATAFOLDER/pileups/$PREFIX.mpileup

