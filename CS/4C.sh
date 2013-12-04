#$ -S /bin/bash
#$ -l h_vmem=10G
#$ -l tmem=10G
#$ -l h_rt=1000:0:0
#$ -cwd
#$ -j y

#tools
TRIMFQ="/home/rmgzshd/symlinks/cohesin-bio/SCRIPT_SV/4C_01_INITIALIZATION_02_trimFASTQ.pl"

#config files
4CSEQPIPECONF="/home/rmgzshd/Sevil4C/4cseqpipe.conf"

export TRIMFQ;

gzip –cd 1_S1_L001_R1_001.fastq.gz  > 1_S1_L001_R1_001.fastq
$TRIMFQ 1_S1_L001_R1_001.fastq
