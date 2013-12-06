#$ -S /bin/bash
#$ -l h_vmem=10G
#$ -l tmem=10G
#$ -l h_rt=1000:0:0
#$ -cwd
#$ -j y

#data directory
DATADIR="/home/rmgzshd/Sevil4C"

# tools
SEQPIPE4C="/SAN/biomed/biomed13/cohesin-bio/4C/4Cseqpipe/src/4cseqpipe.pl"

export DATADIR; export SEQPIPE4C; 

for i in {1:15}
do
	perl $SEQPIPE4C -fastq2raw -ids $i -fastq_fn $DATADIR/1_S1_L001_R1_001.trim.fastq
done
