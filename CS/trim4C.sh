

#tools
TRIMFQ="/home/rmgzshd/symlinks/cohesin-bio/SCRIPT_SV/4C_01_INITIALIZATION_02_trimFASTQ.pl"


export TRIMFQ;

gunzip $DATADIR/1_S1_L001_R1_001.fastq.gz
$TRIMFQ $DATADIR/1_S1_L001_R1_001.fastq
