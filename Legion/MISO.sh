#!/bin/bash

set -o nounset
set -o errexit
set -o pipefail

BAMDATA="/home/rmgzshd/Scratch/Pablo/BAM"
INDEXPY="/shared/ucl/apps/miso/0.5.2/lib/python2.7/site-packages/misopy/index_gff.py"

module unload compilers/intel/11.1/072
module load compilers/gnu/4.6.3
module load python/enthought/7.3-2_2013-10-04
module load samtools/0.1.19
module load bedtools/2.17.0 
module load miso/0.5.2
export ENSEMBL75; export HG19GFF3; export BAMDATA;

python run_events_analysis.py --compute-genes-psi ./indexed $BAMDATA/my_sample1.bam \
--output-dir output \
--read-len 100 \
--use-cluster \

