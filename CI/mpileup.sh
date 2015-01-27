#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

#shortcuts to tools
SAMTOOLS="/home/rmgzshd/samtools/samtools"

# shortcuts to annotation
GENOME="/mnt/store1/GATK_BUNDLE_2.8_hg19/ucsc.hg19.fasta"
CAPTURE_REGION="/mnt/store1/WXS_CAPTURE_BEDS/SeqCap_EZ_Exome_v3_primary.bed"
DATAFOLDER="/mnt/store1/CESC_WXS"

hostname
date
$SAMTOOLS mpileup -A -B -C 50 -d 10000 -m 3 -F 0.0002 -q 20 -Q 20 \
-f $GENOME \
-l $CAPTURE_REGION \
$DATAFOLDER/191e1f11c562128d1f2afb31f50e73b5.bam \
$DATAFOLDER/2308c254ed02bbab73f9478cc398e42b.bam > \
$DATAFOLDER/pileups/TCGA-2W-A8YY.mpileup

