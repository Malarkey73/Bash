#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

# This is all adapted from a Stephen Turner blogpost:
# http://www.gettinggeneticsdone.com/2014/03/visualize-coverage-exome-targeted-ngs-bedtools.html
# I have added a script

#tools
BEDTOOLS="/home/rmgzshd/bedtools2/bin/bedtools"

#places
CAPTURE_REGION="/mnt/store1/WXS_CAPTURE_BEDS/SeqCap_EZ_Exome_v2.target.bed"

# Set up to run in GNU parallel.
find *bam | parallel "$BEDTOOLS coverage -hist -abam {} -b $CAPTURE_REGION | grep ^all > {}.hist.all.txt"

mkdir -p COVERAGE 
mv *.hist.all.txt COVERAGE


# NB There is an RScript on the same page that will view the results
# ANd I have copied that R Script on my github RScipt page

