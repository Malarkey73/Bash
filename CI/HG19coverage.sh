#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail


SAMBAMBA="/home/rmgzshd/sambamba/sambamba"
BEDTOOLS="/home/rmgzshd/bedtools2/bin/bedtools"
LOGFILE="/mnt/store1/cghub/HNSC/coverage.log"
GENOMESIZES="/mnt/store1/cghub/HNSC/Tumour/hg19.chrom.sizes"

REGIONS="2:100000-1000000,3:100000:1000000,4:100000:1000000"


FOUND=$(find -maxdepth 2 -name *.sorted.bam)

for BAM in $FOUND
do
	# calculate coverage NB NEED TO WRITE/APPEND THIS TO A LOG FILE
	PREFIX=$(echo ${BAM} | sed 's/.sorted.bam//')

    HG19COV=$(
    	$SAMBAMBA view --format=bam -F  $BAM $REGIONS 	|
    	$BEDTOOLS genomecov -bg -ibam stdin -g $GENOME  |
   		awk '{sum+=  qq$4} END { print sum/NR}'
   		)

	printf "HG19\t$PREFIX\t$HG19COV\n" >> $LOGFILE

done
