#!/bin/bash

VARSCAN="/home/rmgzshd/VarScan/VarScan.v2.3.7.jar"
GENOME="/mnt/store1/GATK_BUNDLE_2.8_hg19/Homo_sapiens_assembly19.fasta"
BAMDATA="/mnt/store1/CESC_WXS/TCGA-C5-1AMF.bam"

INDELDATA="/mnt/store1/CESC_WXS/TCGA-C5-1AMF.indel"

export VARSCAN; export MPILEUPDATA; export INDELDATA;

/home/gswgoh/software/bam-readcount/bin/bam-readcount $BAMDATA -q 20 -b 20 -f$GENOME -l ~/scratch/varscan/IBC-007T1.snp.Somatic.pos -w 1 > ~/scratch/varscan/IBC-007T1.snp.Somatic.readCount; 
/usr/bin/perl /home/gswgoh/software/VarScan2_scripts/fpfilter.pl ~/scratch/varscan/IBC-007T1.snp.Somatic.pos ~/scratch/varscan/IBC-007T1.snp.Somatic.readCount --output-basename ~/scratch/varscan/IBC-007T1.snp.Somatic