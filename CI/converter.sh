#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

# usage: nohup ./converter.sh S1.vcf & 
# alt usage: ls *.vcf | parallel ./converter.sh {} & 

CONVERT2ANNOVAR="/home/rmgzshd/annovar/convert2annovar.pl"
FILE=`echo $1`
PREFIX=`echo "${FILE%%.*}"`

export CONVERT2ANNOVAR; export FILE; export PREFIX

$CONVERT2ANNOVAR -format vcf4 $FILE > $PREFIX.avinput
