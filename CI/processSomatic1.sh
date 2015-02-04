#!/bin/bash

VARSCAN="/home/rmgzshd/VarScan/VarScan.v2.3.7.jar"
SNPDATA="/mnt/store1/CESC_WXS/TCGA-C5-1AMF.snp"
INDELDATA="/mnt/store1/CESC_WXS/TCGA-C5-1AMF.indel"

export VARSCAN; export MPILEUPDATA; export INDELDATA;

java -d64 -Xmx4g -jar $VARSCAN processSomatic $SNPDATA 

grep Somatic  "$SNPDATA.Somatic" | grep -v _ | awk '{print $1, "\t", $2, "\t", $2}' >  "$SNPDATA.Somatic.pos" 

java -d64 -Xmx4g -jar $VARSCAN processSomatic $INDELDATA
