#!/bin/bash
PREFIX=`echo $1`

cat *_R1_* | zcat -c | gzip -f > $PREFIX.R1.fastq.gz
cat *_R2_* | zcat -c | gzip -f > $PREFIX.R2.fastq.gz