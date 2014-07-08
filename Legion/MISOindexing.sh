#!/bin/bash -l
#$ -S /bin/bash
# 2. Request ten hours of wallclock time (format hours:minutes:seconds).
#$ -l h_rt=2:00:0
# 3. Request 1 gigabyte of RAM 
#$ -l mem=32G
# 4. Request 15 gigabyte of TMPDIR space (default is 10 GB)
#$ -l tmpfs=32G
# 5. Set the name of the job.
#$ -N MISOindexing
# 6. Find <your_project_id> by running the command "groups"
#$ -P TCGAHPVHNSC
# 7. Replace "<your_UCL_id>" with your UCL user ID :)
#$ -wd /home/rmgzshd/Scratch/Pablo/MISO
# 8. request 16 threads
#$ -l thr=16

set -o nounset
set -o errexit
set -o pipefail

HG19GFF3="/home/rmgzshd/Scratch/Genome/MISO_hg19_altevents"
INDEXPY="/shared/ucl/apps/miso/0.5.2/lib/python2.7/site-packages/misopy/index_gff.py"

module unload compilers/intel/11.1/072
module load compilers/gnu/4.6.3
module load python/enthought/7.3-2_2013-10-04
module load samtools/0.1.19
module load bedtools/2.17.0 
module load miso/0.5.2
export ENSEMBL75; export HG19GFF3; export INDEXPY;

python $INDEXPY --index $HG19GFF3/ALL.hg19.gff3 indexed/
