#!/bin/bash

VARSCAN="/home/rmgzshd/VarScan/VarScan.v2.3.7.jar"
MPILEUPDATA="/mnt/store1/CESC_WXS/pileups/TCGA-C5-1AMF.mpileup"
BASENAME="TCGA-C5-1AMF"

export VARSCAN; export MPILEUPDATA; export BASENAME;


java -d64 -Xmx4g -jar $VARSCAN somatic $MPILEUPDATA $BASENAME \
 --mpileup 1 --min-coverage 8 --min-coverage-normal 10 \
 --min-coverage-tumor 6 --min-var-freq 0.01 \
 --min-freq-for-hom 0.75 --normal-purity 1 --p-value 0.99 \
 --somatic-p-value 0.05 --tumor-purity 0.5 --strand-filter 0 