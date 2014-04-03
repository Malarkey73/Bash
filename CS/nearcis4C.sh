#$ -S /bin/bash
#$ -l h_vmem=10G
#$ -l tmem=10G
#$ -l h_rt=1000:0:0
#$ -cwd
#$ -j y

# tools
SEQPIPE4C="/SAN/biomed/biomed13/cohesin-bio/SCRIPTS/4C/4Cseqpipe/4cseqpipe.pl"

export SEQPIPE4C;

perl $SEQPIPE4C -nearcis -calc_from 72500000 -calc_to 73500000 -stat_type median -trend_resolution 5000 -ids 1 -figure_fn AstEsco_raw_Aug13_1.png

perl $SEQPIPE4C -nearcis -calc_from 72500000 -calc_to 73500000 -stat_type median -trend_resolution 5000 -ids 2 -figure_fn AstEsco_raw_Aug13_2.png

