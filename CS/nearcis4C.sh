#$ -S /bin/bash
#$ -l h_vmem=10G
#$ -l tmem=10G
#$ -l h_rt=1000:0:0
#$ -cwd
#$ -j y


# tools
SEQPIPE4C="/home/rmgzshd/FourC/4cseqpipe.pl"

export DATADIR; export SEQPIPE4C; 

perl $SEQPIPE4C -nearcis -calc_from 74000000 -calc_to 76000000 -stat_type median -trend_resolution 5000 -ids 1 -figure_fn NS_2276.png

perl $SEQPIPE4C -nearcis -calc_from 74000000 -calc_to 76000000 -stat_type median -trend_resolution 5000 -ids 2 -figure_fn NS_2282.png

perl $SEQPIPE4C -nearcis -calc_from 74500000 -calc_to 76500000 -stat_type median -trend_resolution 5000 -ids 3 -figure_fn NS_2286.png

perl $SEQPIPE4C -nearcis -calc_from 112000000 -calc_to 114000000 -stat_type median -trend_resolution 5000 -ids 4 -figure_fn NS_14609.png

perl $SEQPIPE4C -nearcis -calc_from 112000000 -calc_to 114000000 -stat_type median -trend_resolution 5000 -ids 5 -figure_fn NS_14610.png

perl $SEQPIPE4C -nearcis -calc_from 112000000 -calc_to 114000000 -stat_type median -trend_resolution 5000 -ids 6 -figure_fn NS_14612.png

perl $SEQPIPE4C -nearcis -calc_from 112000000 -calc_to 114000000 -stat_type median -trend_resolution 5000 -ids 7 -figure_fn NS_14613.png

perl $SEQPIPE4C -nearcis -calc_from 91000000 -calc_to 93000000 -stat_type median -trend_resolution 5000 -ids 8 -figure_fn NS_755.png

perl $SEQPIPE4C -nearcis -calc_from 91000000 -calc_to 93000000 -stat_type median -trend_resolution 5000 -ids 9 -figure_fn NS_756.png

perl $SEQPIPE4C -nearcis -calc_from 91000000 -calc_to 93000000 -stat_type median -trend_resolution 5000 -ids 10 -figure_fn NS_763.png

perl $SEQPIPE4C -nearcis -calc_from 91000000 -calc_to 93000000 -stat_type median -trend_resolution 5000 -ids 11 -figure_fn NS_764.png

perl $SEQPIPE4C -nearcis -calc_from 38500000 -calc_to 40500000 -stat_type median -trend_resolution 5000 -ids 12 -figure_fn NS_222.png

perl $SEQPIPE4C -nearcis -calc_from 38500000 -calc_to 40500000 -stat_type median -trend_resolution 5000 -ids 13 -figure_fn NS_224.png

perl $SEQPIPE4C -nearcis -calc_from 38500000 -calc_to 40500000 -stat_type median -trend_resolution 5000 -ids 14 -figure_fn NS_225.png

perl $SEQPIPE4C -nearcis -calc_from 38500000 -calc_to 40500000 -stat_type median -trend_resolution 5000 -ids 15 -figure_fn NS_227.png
