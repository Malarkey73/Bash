#!/bin/bash

# Usage: script file1 file2 genome species
# E.g.  ./Jenner.sh bamf bamfinput hg19 Hs

#shortcuts to tools
UCSCtools="/home/rmgzshd/UCSCtools"
BEDtools="/home/rmgzshd/bedtools2/bin"
GBINS="/home/rmgzshd/bedtools2/genomes/hg19.10bpbins.bed"

# make sure child processes see shortcuts too
export BEDtools; export UCSCtools; export UCSCtools


# script inputs
ChIP=`echo $1`
Input=`echo $2`
Genome=`echo $3`

# chromosome names
chrN="chr1|chr2|chr3|chr4|chr5|chr6|chr7|chr8|chr9|chr10|chr11|chr12|chr13|chr14|chr15|chr16|chr17|chr18|chr19|chr20|chr21|chr22|chrX|chrY"


$UCSCtools/fetchChromSizes $Genome > "$Genome.chrom.sizes"

#ChIP
(count_pos=$($BEDtools/bamToBed -i "$ChIP.bam" |
$BEDtools/slopBed  -i - -g "$Genome.chrom.sizes" -l 113 -r 36 -s |
grep -Ew $chrN  |
tee "$ChIP.bed" |
wc -l           |
awk '{print 1000000/$1}')
$BEDtools/genomeCoverageBed -bg  -i "$ChIP.bed" -scale $count_pos -g "$Genome.chrom.sizes"  |
$BEDtools/mapBed -a $GBINS -b - -c 4 -o mean |
awk '$4 != "."' >  "$ChIP.bedgraph") &
#Input
(count_pos=$($BEDtools/bamToBed -i "$Input.bam" |
$BEDtools/slopBed  -i - -g "$Genome.chrom.sizes" -l 113 -r 36 -s |
grep -Ew $chrN  |
tee "$Input.bed" |
wc -l           |
awk '{print 1000000/$1}')
$BEDtools/genomeCoverageBed -bg  -i "$Input.bed" -scale $count_pos -g "$Genome.chrom.sizes"  |
$BEDtools/mapBed -a $GBINS -b - -c 4 -o mean |
awk '$4 != "."' >  "$Input.bedgraph") &

wait;

$BEDtools/unionBedGraphs -i "$ChIP.bedgraph" "$Input.bedgraph" | awk '{if ($4 > $5) {print $1,$2,$3,$4-$5}}'> "$ChIP$Input.bedgraph"
$UCSCtools/bedGraphToBigWig "$ChIP$Input.bedgraph" "$Genome.chrom.sizes" "$ChIP.bw"
