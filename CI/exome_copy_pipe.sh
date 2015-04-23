#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

#tools
SAMTOOLS="/home/rmgzshd/samtools-1.2/samtools"
VARSCAN="/home/rmgzshd/VarScan/VarScan.v2.3.7.jar"

#places
GENOME="/mnt/store1/GATK_BUNDLE_2.8_hg19/Homo_sapiens_assembly19.fasta"

#input args
TUMORBAM=`echo $1`
REFBAM=`echo $2`
# longest common prefix from the input BAMS
# and remove any non-alphanumeric final chars e.g. TCGA-C5.A1MF- become TCGA-C5.A1MF
PREFIX=`printf "%s\n%s\n" $TUMORBAM $REFBAM | sed -e 'N;s/^\(.*\).*\n\1.*$/\1/' | sed s'/[^a-zA-Z\d\s:]$//'` 
echo $PREFIX

export GENOME; export TUMORBAM; export REFBAM; export PREFIX; export VARSCAN; export SAMTOOLS

hostname
date

echo "creating raw copynumber file"
$SAMTOOLS mpileup -q 1 -f $GENOME $TUMORBAM $REFBAM   | 
awk '{if($4 >= 6) print $0}' | 
awk '{if($7 != 0) print $0}' | 
java -d64 -Xmx8g -jar $VARSCAN copynumber --output-file $PREFIX --mpileup 1

# this last part is not needed for the sequenza pipeline
#echo "creating adjusted preliminary copynumber files \n"
#java -d64 -Xmx8g -jar $VARSCAN copyCaller output.copynumber --output-file $PREFIX.copynumber.called --output-homdel-file $PREFIX.copynumber.called.homdel

mkdir -p SEQUENZA 
mv $PREFIX.copynumber SEQUENZA
echo "done $PREFIX"
