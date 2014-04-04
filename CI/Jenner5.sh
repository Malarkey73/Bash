#!/bin/bash
# Usage: script file1 file2 genome species
# E.g.  ./Jenner.sh bamfilechip bamfileinput hg19 Hs
# The input bam MUST BE SORTED

#shortcuts to tools
UCSCtools="/home/rmgzshd/UCSCtools"
BEDtools="/home/rmgzshd/bedtools2/bin"
MACS="/home/rmgzshd/MACS14/bin/macs14"

# bedtools makewindows -g hg19.chrom.sizes -w 10 > hg19.10bpbins.bed
# this has to be done once for a genome and particular binsize
GBINS="/home/rmgzshd/bedtools2/genomes/hg19.10bpbins.bed"

# make sure child processes see shortcuts too
export BEDtools; export UCSCtools; export UCSCtools; export MACS;

# script inputs
ChIP=`echo $1`
Input=`echo $2`
Genome=`echo $3`
Species=`echo $4`

# chromosome names
chrN="chr1|chr2|chr3|chr4|chr5|chr6|chr7|chr8|chr9|chr10|chr11|chr12|chr13|chr14|chr15|chr16|chr17|chr18|chr19|chr20|chr21|chr22|chrX|chrY"

# start timer for script
STARTTIME=$(date +%s)

# fetch the chromosome size info only if it doesnt exist
if [ ! -f "$Genome.chrom.sizes" ]; then
	$UCSCtools/fetchChromSizes $Genome > "$Genome.chrom.sizes"
fi

#ChIP sorted bam to scaled, chr filtered bed
(
	count_pos=$($BEDtools/bamToBed -i "$ChIP.sorted.bam" 										|
	grep -Ew $chrN																			|
	uniq																						|
	$BEDtools/slopBed  -i - -g "$Genome.chrom.sizes" -l 113 -r 36 -s 							|
	tee "$ChIP.bed" 																			|
	wc -l           																			|
	awk '{print 1000000/$1}')
	$BEDtools/genomeCoverageBed -bg  -i "$ChIP.bed" -scale $count_pos -g "$Genome.chrom.sizes" 	|
	$BEDtools/mapBed -a $GBINS -b - -c 4 -o mean  >  "$ChIP.bedgraph"
)&

#Input sorted bam to scaled, chr filtered bed
(
	count_pos=$($BEDtools/bamToBed -i "$Input.sorted.bam" 										|
	grep -Ew $chrN  																			|
	uniq																						|
	$BEDtools/slopBed  -i - -g "$Genome.chrom.sizes" -l 113 -r 36 -s 							|
	tee "$Input.bed" 																			|
	wc -l           																			|
	awk '{print 1000000/$1}')
	$BEDtools/genomeCoverageBed -bg  -i "$Input.bed" -scale $count_pos -g "$Genome.chrom.sizes" |
	$BEDtools/mapBed -a $GBINS -b - -c 4 -o mean  >  "$Input.bedgraph"
)&

# have to wait for both to finish
wait

#Run two MACS at same time
$MACS -t "$Input.bed"  -c "$ChIP.bed" -g $Species --keep-dup=1 -n "$ChIP"_10-9 -p 1e-6 &
$MACS -t "$Input.bed"  -c "$ChIP.bed" -g $Species --keep-dup=1 -n "$ChIP"_10-7 -p 1e-6 &

# subtract Input from ChIP, where ChIP > Input subtract, 
# otherwise skip that line, then > bedgraph && bigwig
paste "$ChIP.bedgraph" "$Input.bedgraph" | awk '{if ($8 != "." && $4 > $8) {print $1,$2,$3,$4-$8}}'> "$ChIP$Input.bedgraph" &&
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
