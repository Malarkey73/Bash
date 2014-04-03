#$ -S /bin/bash
#$ -l h_vmem=10G
#$ -l tmem=10G
#$ -l h_rt=1000:0:0
#$ -cwd
#$ -j y

#data directory
DATADIR="/SAN/biomed/biomed13/cohesin-bio/SCRIPTS/4C/4Cseqpipe/rawdata"

# tools
SEQPIPE4C="/SAN/biomed/biomed13/cohesin-bio/SCRIPTS/4C/4Cseqpipe/4cseqpipe.pl"

export DATADIR; export SEQPIPE4C; 

perl $SEQPIPE4C -map -ids 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15
