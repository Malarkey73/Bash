#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

BT2="/home/rmgzshd/bowtie2-2.1.0/bowtie2"
BT2GENOME="/mnt/store1/cghub/HPV_BT2/HPV"
SAMBAMBA="/home/rmgzshd/sambamba/sambamba"
BEDTOOLS="/home/rmgzshd/bedtools2/bin/bedtools"
GENOME="/mnt/store1/cghub/HPV_BT2/HPV.sizes"

export BT2; export BT2GENOME; export SAMBAMBA; export BEDTOOLS;

convertsecs() {
 ((h=${1}/3600))
 ((m=(${1}%3600)/60))
 ((s=${1}%60))
 printf "\t%02d:%02d:%02d\n" $h $m $s
}
STARTTIME=$(date +%s)

for BAM in *.sorted.bam
do
	
	PREFIX=$(echo ${BAM} | sed 's/.sorted.bam//')
	printf "\n\n beginning %s sample.\n " $PREFIX
	
	$SAMBAMBA view --format=bam -F  "paired and (unmapped or mate_is_unmapped)" $BAM 				|
	tee >($BEDTOOLS bamtofastq -i stdin -fq R1.fq -fq2 R2.fq)										|
	$BEDTOOLS bamtobed -bedpe -i stdin 	| sort -S8G -k 7											|
	# I clip the mate pair ID last char read flag so I can merge with HPV bed  
	awk '{if ($2 != -1 || $5 != -1) print $1,$2,$3,$4,$5,$6, substr($7, 0, length($7)-2), $8, $9, $10}' > $PREFIX.anchor.bedpe
	
	# print time
	NEXTTIME1=$(date +%s)
#	printf "\n Finished extracting unmapped reads in: \n"
#	convertsecs $(($NEXTTIME1 - $STARTTIME))
		
	# Use bowtie2 to find matches amongst the unmapped genome in HPV 	
	$BT2 -p 24 -x $BT2GENOME -1 R1.fq -2 R2.fq 				| 
	$SAMBAMBA view -S --format=bam /dev/stdin 				|
	tee $PREFIX.hpv.bam 							|
    	$BEDTOOLS bamtobed -bedpe -i stdin | sort -S8G -k 7 			|
	awk '{if ($2 != -1 || $5 != -1) print $1,$2,$3,$4,$5, $6, substr($7, 0, length($7)-2), $8, $9, $10}' > $PREFIX.hpv.bedpe

    	# print time
	NEXTTIME2=$(date +%s)
#	printf "\n Finished aligning to HPV genome in: \n"
#	convertsecs $(($NEXTTIME2 - $NEXTTIME1))
    
    	# calculate coverage NB NEED TO WRITE/APPEND THIS TO A LOG FILE
    	printf "\n HPV coverage :"
    	$BEDTOOLS genomecov -bg -ibam $PREFIX.hpv.bam  -g $GENOME  |
    	awk '{sum+=  qq$4} END { print sum/NR}'
	
	# print time
	NEXTTIME3=$(date +%s)
#	printf "\n Finished calculating HPV coverage in: \n"
#	convertsecs $(($NEXTTIME3 - $NEXTTIME2))

	#awk '{if ($2!= -1 || $5 != -1)  print}'
	#join -1 4 -2 4 sort.hpv.bed sort.anchor.bed > $PREFIX.join.bed

	printf "\n Next file?: "
done
printf "\n No, all done: "

rm R1.fq
rm R2.fq
# print time
ENDTIME=$(date +%s)
#printf "\n\t Finished all work in: "
#convertsecs $(($ENDTIME - $STARTTIME))
