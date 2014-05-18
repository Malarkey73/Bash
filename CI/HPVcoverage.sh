#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

BT2="/home/rmgzshd/bowtie2-2.1.0/bowtie2"
BT2GENOME="/mnt/store1/cghub/HPV_BT2/HPV"
SAMBAMBA="/home/rmgzshd/sambamba/sambamba"
BEDTOOLS="/home/rmgzshd/bedtools2/bin/bedtools"
GENOMESIZES="/mnt/store1/cghub/HPV_BT2/HPV.sizes"

export BT2; export BT2GENOME; export SAMBAMBA; export BEDTOOLS;

convertsecs() {
 ((h=${1}/3600))
 ((m=(${1}%3600)/60))
 ((s=${1}%60))
 printf "%02d:%02d:%02d\n" $h $m $s
}

# start timer for script
STARTTIME=$(date +%s)



for BAM in *.sorted.bam
do
	
	PREFIX=$(echo ${BAM} | sed 's/.sorted.bam//')
	printf "\n\n\t beginning %s \n " $PREFIX
	
	$SAMBAMBA view --format=bam -F  "paired and (unmapped or mate_is_unmapped)" $BAM |
	tee >($BEDTOOLS bamtofastq -i stdin -fq R1.fq -fq2 R2.fq)			|
	$BEDTOOLS bamtobed -i stdin 										|
	# I clip the mate pair char so I can easily compare to HPV bed  
	awk '{print $1,$2,$3,substr($4, 0, length($4)-2),$5, $6}' 			> 
	$PREFIX.anchor.bed
	# real	3m56.793s on 20GB
	# user	17m21.587s
	# sys	0m39.095s

	# print time
	ENDTIME=$(date +%s)
	printf '\n\t Finished extracting unmapped reads in: '
	echo convertsecs $(($ENDTIME - $STARTTIME))
	NEXTTIME==$(date +%s)

	
	# Use bowtie2 to find matches amongst the unmapped genome in HPV 	
	$BT2 -p 24 -x $BT2GENOME -1 R1.fq -2 R2.fq 							| 
	$SAMBAMBA view -S --format=bam /dev/stdin 							|
    $BEDTOOLS bamtobed -i stdin  										|
    # NB same clip as ANCHOR bed above so can compare mate pairs
	awk '{print $1,$2,$3,substr($4, 0, length($4)-2),$5, $6}' 			> 
    $PREFIX.hpv.bed

    # print time
	ENDTIME=$(date +%s)
	printf '\n\t Finished aligning to HPV genome in: '
	echo convertsecs $(($ENDTIME - $NEXTTIME))
	NEXTTIME==$(date +%s)
    
    # calculate coverage
    $BEDTOOLS genomecov -i $PREFIX.hpv.bed -g $GENOME -bg 
	awk '{sum+=$4} END { print sum/NR}'
	# mean depth 138.822 on 20GB HPV16+ HNSC file TCGA-BA-5153

	# print time
	ENDTIME=$(date +%s)
	printf '\n\t Finished calculating HPV coverage in: '
	echo convertsecs $(($ENDTIME - $NEXTTIME))
	NEXTTIME==$(date +%s)


	printf '\n\n\t Next file: '
done

rm R1.fq
rm R2.fq

# print time
ENDTIME=$(date +%s)
printf '\n\t Finished all work in: '
echo convertsecs $(($ENDTIME - $STARTTIME))





