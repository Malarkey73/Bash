#!/bin/bash

set -o nounset
set -o errexit
set -o pipefail

HG19GFF3="/home/rmgzshd/Scratch/Genome/ENSEMBL.homo_sapiens.release-75/hg19"
INDEXPY="/shared/ucl/apps/miso/0.5.2/lib/python2.7/site-packages/misopy/index_gff.py"

module unload compilers/intel/11.1/072
module load compilers/gnu/4.6.3
module load python/enthought/7.3-2_2013-10-04
module load samtools/0.1.19
module load bedtools/2.17.0 
module load miso/0.5.2
export ENSEMBL75; export HG19GFF3; export INDEXPY;

python $INDEXPY --index $HG19GFF3/ALL.hg19.gff3 indexed/
