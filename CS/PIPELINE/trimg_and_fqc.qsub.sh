#!/bin/bash -l
#$ -S /bin/bash
## Time limit
#$ -l h_rt=2:00:00
## Memory requirements (h_vmem and them must match)
#$ -l h_vmem=5G
#$ -l tmem=5G
#$ -l fastio=1
#$ -l scratch0free=10G
## Use the current working directory for running the job
#$ -cwd


# Paths
TRIMGALORE="/share/apps/genomics/trim_galore/trim_galore"



export TRIMGALORE;


# Args
R1PREFIX=$1
PAIRED=$2
GZ=$3
PREFIX=`echo $R1PREFIX | sed -e 's/.R1_001//g'`
R2PREFIX=`echo $R1PREFIX | sed -e 's/R1_001/R2_001/g'`
DATA_DIR=$PWD


# These echo commands are for debugging. If there is a problem
# then check these ar ecorrectly set in the script.oxxxxx
# file produced
echo $R1PREFIX
echo $PAIRED
echo $GZ
echo $PREFIX
echo $R2PREFIX
echo $DATA_DIR


# create new self scrubbing TMPDIR and DATA_DIR/QC for output
mkdir -p /scratch0/$TMPDIR
rmdir $TMPDIR
ln -s /scratch0/$TMPDIR $TMPDIR
mkdir -p $DATA_DIR/FASTQC

# Do this if the fastq are single end
if [ $PAIRED == 0 ]; then
	rsync -aL $DATA_DIR/$R1PREFIX.fastq* $TMPDIR
	cd $TMPDIR
	# if not gzipped
	if [ $GZ == 0 ]; then
		TRIMGALORE --fastqc --dont_gzip --phred33 $R1PREFIX.fastq

	fi
	# if gzipped
	if [ $GZ == 1 ]; then
		TRIMGALORE --fastqc --dont_gzip --phred33 $R1PREFIX.fastq
	fi
	rsync -a $R1PREFIX*.zip $DATA_DIR
	rsync -a $R1PREFIX*.txt $DATA_DIR
	rename _trimmed.fq .fastq $PREFIX_trimmed.fq
	rsync -a $PREFIX*.fastq $DATA_DIR

fi

# Or do this if the fastq are paired. This is not necessary
# here for FastQC but many programs will require dual R1 & R2
# input so this template handles pairs.
# NB note the use of & to bg processes and speed stuff up a bit
# NB note to the use of "wait" to make sure bg commands are done
# before moving on

if [ $PAIRED == 1 ]; then
	rsync -aL $DATA_DIR/$R1PREFIX.fastq* $TMPDIR &
    rsync -aL $DATA_DIR/$R2PREFIX.fastq* $TMPDIR &
    wait
    cd $TMPDIR
    # if not gzipped
	if [ $GZ == 0 ]; then
    	TRIMGALORE --fastqc --dont_gzip --phred33 --paired $R1PREFIX.fastq $R2PREFIX.fastq
	fi
	# if gzipped
	if [ $GZ == 1 ]; then
		TRIMGALORE --fastqc --dont_gzip --phred33 --paired $R1PREFIX.fastq $R2PREFIX.fastq
	fi	

	rsync -a $PREFIX*.zip $DATA_DIR
	rsync -a $PREFIX*.txt $DATA_DIR
	mv $R1PREFIX_val_1.fq $R1PREFIX.fastq
	mv $R2PREFIX_val_2.fq $R2PREFIX.fastq
	rsync -a $PREFIX*.fastq $DATA_DIR
fi

