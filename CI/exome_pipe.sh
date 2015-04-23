#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

#tools
SAMTOOLS="/home/rmgzshd/samtools-1.2/samtools"
VARSCAN="/home/rmgzshd/VarScan/VarScan.v2.3.7.jar"
BAMREADCOUNT="/home/rmgzshd/bam-readcount/bin/bam-readcount"
FPFILTER="/home/rmgzshd/VarScan/fpfilter.pl"

#places
GENOME="/mnt/store1/GATK_BUNDLE_2.8_hg19/Homo_sapiens_assembly19.fasta"
CAPTURE_REGION="/mnt/store1/WXS_CAPTURE_BEDS/SeqCap_EZ_Exome_v2.target.bed"

#input args
TUMORBAM=`echo $1`
REFBAM=`echo $2`
# longest common prefix from the input BAMS
# and remove any non-alphanumeric final chars e.g. TCGA-C5.A1MF- become TCGA-C5.A1MF
PREFIX=`printf "%s\n%s\n" $TUMORBAM $REFBAM | sed -e 'N;s/^\(.*\).*\n\1.*$/\1/' | sed s'/[^a-zA-Z\d\s:]$//'` 
echo $PREFIX

export GENOME; export CAPTURE_REGION; export DATAFOLDER; export TUMORBAM; 
export REFBAM; export PREFIX; export VARSCAN; export SAMTOOLS

hostname
date
# I changed this to $PREFIX.fifo so that I could run this script in parallel with another
# i.e. there would be 2 separate named fifo
mkfifo $PREFIX.fifo

# th epileup is in a fifo so there is no intermediate file
echo "piling up. \n"
$SAMTOOLS mpileup -A -B -C 50 -d 10000 -m 3 -F 0.0002 -q 20 -Q 20 \
-f $GENOME \
-l $CAPTURE_REGION \
$REFBAM $TUMORBAM > $PREFIX.fifo &

# the inetrmediary pipe is fed to Varscan Sopmatic 
java -d64 -Xmx8g -jar $VARSCAN somatic $PREFIX.fifo $PREFIX \
--mpileup 1 --min-coverage 8 --min-coverage-normal 10 --min-coverage-tumor 6 \
--min-var-freq 0.01 --min-freq-for-hom 0.75 --normal-purity 1 --p-value 0.99 \
--somatic-p-value 0.05 --tumor-purity 0.5 --strand-filter 0

wait
echo "processing step 1. \n"
rm $PREFIX.fifo
java -d64 -Xmx8g -jar $VARSCAN processSomatic $PREFIX.snp &
java -d64 -Xmx8qg -jar $VARSCAN processSomatic $PREFIX.indel &

wait
echo "processing step 2. \n"

grep Somatic $PREFIX.snp.Somatic | grep -v _ | awk '{print $1, "\t", $2, "\t", $2}' >  $PREFIX.snp.Somatic.pos
$BAMREADCOUNT $TUMORBAM -q 20 -b 20 -f $GENOME -l $PREFIX.snp.Somatic.pos -w 1 > $PREFIX.snp.Somatic.rc 
/usr/bin/perl $FPFILTER $PREFIX.snp.Somatic $PREFIX.snp.Somatic.rc --output-basename $PREFIX.snp.Somatic

# this makes the input format ANNOVAR like
awk  '{print $1,$2,$2,$3,$4}' $PREFIX.snp.Somatic.pass > temp
temp > $PREFIX.snp.Somatic.pass

rm temp
rm $PREFIX.snp.Somatic.pos
rm $PREFIX.snp.Somatic.rc
rm $PREFIX.snp.Somatic
rm $PREFIX.snp.hc
rm $PREFIX.snp.Germline.hc
rm $PREFIX.snp.LOH*


echo "done"
