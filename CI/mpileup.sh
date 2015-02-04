#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

#shortcuts to tools
SAMTOOLS="/home/rmgzshd/samtools/samtools"

# shortcuts to annotation
GENOME="/mnt/store1/GATK_BUNDLE_2.8_hg19/Homo_sapiens_assembly19.fasta"
CAPTURE_REGION="/mnt/store1/WXS_CAPTURE_BEDS/whole_exome_agilent_1.1_refseq_plus_3_boosters.targetIntervals.bed"
DATAFOLDER="/mnt/store1/CESC_WXS"
TUMORBAM=`echo $1`
REFBAM=`echo $2`
# longest common prefix
PREFIX=`printf "%s\n%s\n" $TUMORBAM $REFBAM | sed -e 'N;s/^\(.*\).*\n\1.*$/\1/'`
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

