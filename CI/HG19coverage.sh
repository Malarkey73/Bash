#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail


SAMBAMBA="/home/rmgzshd/sambamba/sambamba"
BEDTOOLS="/home/rmgzshd/bedtools2/bin/bedtools"
LOGFILE="/mnt/store1/cghub/HNSC/Tumour/coverage.log"
GENOMESIZES="/mnt/store1/cghub/HNSC/Tumour/hg19.chrom.sizes"
REGIONS="2:100000-1000000,3:100000-1000000,4:100000-1000000,5:100000-1000000"

export SAMBAMBA; export BEDTOOLS; export LOGFILE; export GENOMESIZES; export REGIONS;

# start timer for script
STARTTIME=$(date +%s)

FOUND=$(find -maxdepth 2 -name *.sorted.bam)

for BAM in $FOUND
do
	# calculate coverage NB NEED TO WRITE/APPEND THIS TO A LOG FILE
	PREFIX=$(echo ${BAM} | sed 's/.sorted.bam//')

    HG19COV=$(
    	$SAMBAMBA view --format=bam $BAM $REGIONS 		|
    	$BEDTOOLS genomecov -bg -ibam stdin -g $GENOMESIZES  |
   		awk '{sum+=  $4} END { print sum/NR}'
   		)

	printf "HG19\t$PREFIX\t$HG19COV\n" >> $LOGFILE

done

# print time, tested ca 8 mins for 2 * 20GB files
ENDTIME=$(date +%s)
dt=$(($ENDTIME - $STARTTIME))
ds=$((dt % 60))
dm=$(((dt / 60) % 60))
dh=$((dt / 3600))
printf '\n\t %d:%02d:%02d \n' $dh $dm $ds

