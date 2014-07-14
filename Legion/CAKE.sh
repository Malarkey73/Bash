#!/bin/bash -l
#$ -S /bin/bash
# 2. Request ten hours of wallclock time (format hours:minutes:seconds).
#$ -l h_rt=10:00:0
# 3. Request 1 gigabyte of RAM 
#$ -l mem=32G
# 4. Request 15 gigabyte of TMPDIR space (default is 10 GB)
#$ -l tmpfs=32G
# 5. Set the name of the job.
#$ -N CAKEHPV
# 6. Find <your_project_id> by running the command "groups"
#$ -P TCGAHPVHNSC
# 7. Replace "<your_UCL_id>" with your UCL user ID :)
#$ -wd /home/rmgzshd/Scratch/CAKE/BAM
# 8. request 16 threads
#$ -l thr=16

set -o nounset
set -o errexit
set -o pipefail

RUNCAKE="/shared/ucl/apps/cake/1.0/trunk/scripts/run_somatic_pipeline.pl"
module unload compilers/intel/11.1/072
module load compilers/gnu/4.6.3
module load samtools/0.1.19
module load cake/1.0

export ENSEMBL75; export FQDATA; export SEQTK;

perl $RUNCAKE \
     -s /home/rmgzshd/Scratch/CAKE/CAKE.sif \
     -species human \
     -callers mpileup,varscan,bambino,somaticsniper \
     -separator \t \
     -o /home/rmgzshd/Scratch/CAKE \
     -mode [FULL|CALLING|FILTERING]

