#!/bin/bash
# Usage: script file1 file2 genome species
# E.g.  ./Jenner.sh bamfilechip bamfileinput hg19 Hs
# The input bam MUST BE SORTED
# NB the genome chrom.sizes if it already exist MUST BE SORTED IN THE SAME ORDER e.g. chr1:22,X,Y

#shortcuts to tools
UCSCtools="/home/rmgzshd/UCSCtools"
BEDtools="/home/rmgzshd/bedtools2/bin"
MACS="/home/rmgzshd/MACS14/bin/macs14"

# NB This MUST BE CREATED from a genome chr1:22,X,Y too
# e.g. bedtools makewindows -g hg19.chrom.sizes -w 10 > hg19.10bpbins.bed
GBINS="/home/rmgzshd/bedtools2/genomes/hg19.10bpbins.bed"

# make sure child processes see shortcuts too
export BEDtools; export UCSCtools; export UCSCtools; export MACS;

# script input arguments
ChIP=`echo $1`
Input=`echo $2`
Genome=`echo $3`
Species=`echo $4`

# chromosome names
chrN="chr1|chr2|chr3|chr4|chr5|chr6|chr7|chr8|chr9|chr10|chr11|chr12|chr13|chr14|chr15|chr16|chr17|chr18|chr19|chr20|chr21|chr22|chrX|chrY"
# start timer for script
STARTTIME=$(date +%s)

# UCSC fetches chromosomes sorted by size bur we need them chr1:22,X,Y like samtools sort
if [ ! -f "$Genome.chrom.sizes" ]; then
	$UCSCtools/fetchChromSizes hg19 |grep -Ew $chrN | sed 's/chr/chr /' | 
	sed 's/ X/ 23/' | sed 's/ Y/ 24/' | sort -k2,2n | sed 's/chr 23/chrX/' | 
	sed 's/chr 24/chrY/' | sed 's/chr /chr/' > "$Genome.chrom.sizes"
fi

# briefly line by line: bam to bed then
# data filtered down to main chromosome
# only unique reads
# adjust the ends
# tee off a bed file
# whilst counting the filtered unique reads
# save that as rpm scaling factor count_pos
# calculate coverage and pipe into
# 10 bp bins (pre calculated reference)
(
	count_pos=$($BEDtools/bamToBed -i "$ChIP.sorted.bam" |
	grep -Ew $chrN  |
	awk '!array[$2,$3,$6]++' |
	$BEDtools/slopBed  -i - -g "$Genome.chrom.sizes" -l 113 -r 36 -s |
	tee "$ChIP.bed" |
	wc -l |
	awk '{print 1000000/$1}')
	$BEDtools/genomeCoverageBed -bg  -i "$ChIP.bed" -scale $count_pos -g "$Genome.chrom.sizes" |
	$BEDtools/mapBed -a $GBINS -b - -c 4 -o mean -g "$Genome.chrom.sizes" >  "$ChIP.bedgraph"
)&

#same for input consecutively
(
	count_pos=$($BEDtools/bamToBed -i "$Input.sorted.bam" |
	grep -Ew $chrN 	|
	awk '!array[$2,$3,$6]++' |
	$BEDtools/slopBed  -i - -g "$Genome.chrom.sizes" -l 113 -r 36 -s |
	tee "$Input.bed" |
	wc -l |
	awk '{print 1000000/$1}')
	$BEDtools/genomeCoverageBed -bg  -i "$Input.bed" -scale $count_pos -g "$Genome.chrom.sizes" |
	$BEDtools/mapBed -a $GBINS -b - -c 4 -o mean -g "$Genome.chrom.sizes" >  "$Input.bedgraph"
)&

# have to wait for both to finish
wait

#Run two MACS at same time
#$MACS -t "$Input.bed"  -c "$ChIP.bed" -g $Species --keep-dup=1 -n "$ChIP"_10-9 -p 1e-6 &
#$MACS -t "$Input.bed"  -c "$ChIP.bed" -g $Species --keep-dup=1 -n "$ChIP"_10-7 -p 1e-6 &

# Where ChIP > Input , ChIP - Input, otherwise 0
# then to bigwig
paste "$ChIP.bedgraph" "$Input.bedgraph" | awk '{if ($8 != "." && $4 > $8) {print $1,$2,$3,$4-$8}else{print $1,$2,$3,0}}' > "$ChIP$Input.bedgraph" &&
$UCSCtools/bedGraphToBigWig "$ChIP$Input.bedgraph" "$Genome.chrom.sizes" "$ChIP.bw"

wait

# tidy up
rm "$ChIP.bed"
rm "$Input.bed"
rm "$ChIP$Input.bedgraph"
rm "$ChIP.bedgraph"
rm "$Input.bedgraph"

# print time
ENDTIME=$(date +%s)
dt=$(($ENDTIME - $STARTTIME))
ds=$((dt % 60))
dm=$(((dt / 60) % 60))
dh=$((dt / 3600))
printf '\n\t %d:%02d:%02d \n' $dh $dm $ds
