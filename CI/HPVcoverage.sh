#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

# profiled 15GB file in:
#real	11m23.279s
#user	70m25.134s
#sys	9m2.536s

#example log output
# HPV	TCGA-BA-5153-01A-01D-1431_120423_SN590_0154_AC0JBHACXX_s_5_rg	138.75

BT2="/home/rmgzshd/bowtie2-2.1.0/bowtie2"
BT2GENOME="/mnt/store1/cghub/HPV_BT2/HPV"
SAMBAMBA="/home/rmgzshd/sambamba/sambamba"
SAMTOOLS="/home/rmgzshd/samtools/samtools"
BEDTOOLS="/home/rmgzshd/bedtools2/bin/bedtools"
GENOMESIZES="/mnt/store1/cghub/HPV_BT2/HPV.sizes"
LOGFILE="/mnt/store1/cghub/HNSC/Tumour/coverage.log"
RESULTS="/mnt/store1/cghub/HNSC/Tumour/Results"

export BT2; export BT2GENOME; export SAMBAMBA; export BEDTOOLS; export RESULTS; export BAMTOOLS;

convertsecs() {
 ((h=${1}/3600))
 ((m=(${1}%3600)/60))
 ((s=${1}%60))
 printf "\t%02d:%02d:%02d\n" $h $m $s
}
STARTTIME=$(date +%s)


FOUND=$(find -maxdepth 2 -name *.sorted.bam)

for BAM in $FOUND
do
	
	PREFIX=$(echo ${BAM} | sed 's/.sorted.bam//')
	printf "\n\n beginning %s sample.\n " $PREFIX
	
	$SAMBAMBA view --format=bam -F  "paired and (unmapped or mate_is_unmapped)" $BAM 				|
	tee >($BEDTOOLS bamtofastq -i stdin -fq R1.fq -fq2 R2.fq)										|
	$BEDTOOLS bamtobed -bedpe -i stdin 	| sort -S8G -k 7											|
	# I clip the mate pair ID last char read flag so I can merge with HPV bed  
	#awk '{if ($2 != -1 || $5 != -1) print $1,$2,$3,$4,$5,$6, substr($7, 0, length($7)-2), $8, $9, $10}' > $PREFIX.anchor.bedpe
	awk '{if ($2 != -1 || $5 != -1) print $1,$2,$3,$4,$5,$6, $7, $8, $9, $10}' > $PREFIX.anchor.bedpe
	

	# print time
	NEXTTIME1=$(date +%s)
	printf "\n Finished extracting unmapped reads: \n"
	#	convertsecs $(($NEXTTIME1 - $STARTTIME))
		
	# Use bowtie2 to find matches amongst the unmapped genome in HPV
	# convert to bam
	# sort it 	
	$BT2 -p 24 -x $BT2GENOME -1 R1.fq -2 R2.fq 				|
	$SAMBAMBA view -S --format=bam /dev/stdin > tempsortbam	 
	#I have had problems with stability of sambabmba sort of BIG files (>50GB)
	#$SAMBAMBA sort -m=32G tempsortbam -o $PREFIX.hpv.bam 
	$SAMTOOLS sort -@ 12 -m 32G tempsortbam $PREFIX.hpv
	
	# NB sambabmba sort whilst fast doesn't work on stream???
	# So it fucks up the pipe flow here (see above)
	$BEDTOOLS bamtobed -bedpe -i $PREFIX.hpv.bam | sort -S8G -k 7 			|
	#awk '{if ($2 != -1 || $5 != -1) print $1,$2,$3,$4,$5, $6, substr($7, 0, length($7)-2), $8, $9, $10}' > $PREFIX.hpv.bedpe		awk '{if ($2 != -1 || $5 != -1) print $1,$2,$3,$4,$5, $6, substr($7, 0, length($7)-2), $8, $9, $10}' > $PREFIX.hpv.bedpe
	awk '{if ($2 != -1 || $5 != -1) print $1,$2,$3,$4,$5, $6, $7, $8, $9, $10}' > $PREFIX.hpv.bedpe

   	# print time
	NEXTTIME2=$(date +%s)
	printf "\n Finished aligning to HPV genome: \n"
	#	convertsecs $(($NEXTTIME2 - $NEXTTIME1))
    
    # calculate coverage NB NEED TO WRITE/APPEND THIS TO A LOG FILE
    printf "\n HPV coverage :"
    HPVCOV=$($BEDTOOLS genomecov -bg -ibam $PREFIX.hpv.bam  -g $GENOMESIZES  |
   	awk '{sum+=  $4} END { print sum/NR}')
	printf "HPV\t$PREFIX\t$HPVCOV\n" >> $LOGFILE

	# print time
	NEXTTIME3=$(date +%s)
	printf "\n Finished calculating HPV coverage: \n"
	#	convertsecs $(($NEXTTIME3 - $NEXTTIME2))

	# merge the single mapped HPV and ANCHOR samples to look for bridges
	join -1 7 -2 7 $PREFIX.hpv.bedpe $PREFIX.anchor.bedpe | 
	sort -k1,14 -k2,15 >  $PREFIX.integration.bed

	# print time
	NEXTTIME4=$(date +%s)
	printf "\n Finished integration pairs join: \n"
	#	convertsecs $(($NEXTTIME4 - $NEXTTIME3))

	rm $PREFIX.hpv.bedpe 
	rm $PREFIX.anchor.bedpe
	mv $PREFIX.hpv.bam $RESULTS
	mv $PREFIX.integration.bed $RESULTS
	printf "\n Next file?: "
done
printf "\n No, all done! \n"


rm R1.fq
rm R2.fq
rm tempsortbam

# print time
ENDTIME=$(date +%s)
#printf "\n\t Finished all work in: "
#convertsecs $(($ENDTIME - $STARTTIME))
