#!/bin/bash

# test parameters args
PROGRAM=$1
if [ -n "$PROGRAM" ]; then
        echo "Program is $PROGRAM"
fi




# BAM or FASTQ
FQn=`ls *.fastq* | wc -l`
BAMn=`ls *.bam | wc -l`

# If it is fastq files
if [ $FQn > 0 ]; then
	echo "Total Fastq files= $FQn"

	# detect correct pairing, 0 unpaired, 1 paired
	PAIRED=0
	PAIR1=`ls *R1_001.fastq* | wc -l`
	PAIR2=`ls *R2_001.fastq* | wc -l`
	echo "Total R1s = $PAIR1"
	echo "Total R2s = $PAIR2"
	if [ "$PAIR1" > 0 ] && [ "$PAIR1" == "$PAIR2" ]; then
        echo "R1s and R2s correctly paired."
        PAIRED=1
	fi

	GZ=0
	for FQ in *R1_001.fastq*
	do
		if [[ $FQ =~ \.gz$ ]]; then
			GZ=1
			echo "Data is gzipped."
		fi
		# the prefix of those fastq files
		PREFIX=`echo $FQ | sed -e 's/.fastq$//g;s/.fastq.gz$//g'`
		echo "Prefix is $PREFIX"
		echo "qsub command is qsub $PROGRAM $PREFIX $PAIRED $GZ"
		qsub $PROGRAM $PREFIX $PAIRED $GZ

	done
fi

# If it is bam files (its simpler)
if [ $BAMn > 0 ]; then
	for BAM in *.bam
	do
		PREFIX=`echo $BAM | sed -e 's/.bam$//g'`
		echo "Prefix is $PREFIX"
		echo "qsub command is qsub $PROGRAM $PREFIX "
		qsub $PROGRAM $PREFIX
	done
fi

