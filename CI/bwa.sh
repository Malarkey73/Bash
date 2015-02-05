#!/bin/bash

GENOME=""
PE1=`echo $1`
PE2=`echo $2`
PREFIX=`printf "%s\n%s\n" $PE1 $PE2 | sed -e 'N;s/^\(.*\).*\n\1.*$/\1/' | sed s'/[^a-zA-Z\d\s:]$//'` 

bwa mem -t 4 -M -O

bwa sampe

java -Xmx2g -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R <ref.fa> -I <lane.bam> -o <lane.intervals> --known <bundle/b38/Mills1000G.b38.vcf>
java -Xmx4g -jar GenomeAnalysisTK.jar -T IndelRealigner -R <ref.fa> -I <lane.bam> -targetIntervals <lane.intervals> --known <bundle/b38/Mills1000G.b38.vcf> -o <lane_realigned.bam>