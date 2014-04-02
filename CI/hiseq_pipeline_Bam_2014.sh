#!/bin/bash

#===================================================================================
# ChIP-seq pipeline from sorted.Bam - Richard Jenner 24-1-14
# This pipeline will create BigWig files 
# and run MACS at 2 different p-value cutoffs.
# Uses Jan 2014 update for Hg19 refFlat
# Gave tmp files more useful names
#===================================================================================
lastname=`echo $1`
lastname1=`echo $2`
#===================================================================================
# Converting bam to bedgraph
#===================================================================================
# Treatment file
#===================================================================================
fetchChromSizes $3 > $3".chrom.sizes" 

bamToBed -i $lastname".sorted.bam" > $lastname".bed"

perl bed2solexa_result.pl $lastname".bed" > $lastname".solexa_result"
perl uniq_tags.pl $lastname".bed" > $lastname"_uniq.bed"
bedSort $lastname"_uniq.bed" $lastname"_uniq.bed"
rm tmp"_uniq.bed"
touch tmp"_uniq.bed"
for file1 in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY
do
	grep -w $file1 $lastname"_uniq.bed" >> tmp"_uniq.bed"
done

mv tmp"_uniq.bed" $lastname"_uniq.bed"
count_pos=`wc -l $lastname"_uniq.bed"|awk '{print 1000000/$1}'`
genomeCoverageBed -bg  -i $lastname"_uniq.bed" -scale $count_pos -g $3".chrom.sizes"  > ChIP_tmp1.bedgraph
perl bedgraph_10_windows.pl $3".chrom.sizes" ChIP_tmp1.bedgraph > ChIP_tmp2.bedgraph
#===================================================================================
# Control file
#===================================================================================
bamToBed -i $lastname1".sorted.bam" > $lastname1".bed"
perl bed2solexa_result.pl $lastname".bed" > $lastname".solexa_result"
perl uniq_tags.pl $lastname1".bed" > $lastname1"_uniq.bed"
bedSort $lastname1"_uniq.bed" $lastname1"_uniq.bed"
rm tmp"_uniq.bed"
touch tmp"_uniq.bed"
for file1 in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY
do
	grep -w $file1 $lastname1"_uniq.bed" >> tmp"_uniq.bed"
done
mv tmp"_uniq.bed" $lastname1"_uniq.bed"
count_pos=`wc -l $lastname1"_uniq.bed"|awk '{print 1000000/$1}'`
genomeCoverageBed -bg  -i $lastname1"_uniq.bed" -scale $count_pos -g $3".chrom.sizes"  > Input_tmp1.bedgraph
perl bedgraph_10_windows.pl $3".chrom.sizes" Input_tmp1.bedgraph > Input_tmp2.bedgraph
bedGraphToBigWig Input_tmp2.bedgraph $3".chrom.sizes" $lastname1".bw"
#===================================================================================
# Background corrected Bigwig
#===================================================================================
unionBedGraphs -i ChIP_tmp2.bedgraph Input_tmp2.bedgraph | awk '{if ($4 > $5) {print $1,$2,$3,$4-$5}else{print $1,$2,$3,0}}'| sed -e 's/ /\t/g' > ChIP-Input_tmp.bedgraph
#perl bedgraph_10_windows.pl $3".chrom.sizes" tmp2.bedgraph > tmp21.bedgraph
bedGraphToBigWig ChIP-Input_tmp.bedgraph $3".chrom.sizes" $lastname".bw"

#===================================================================================
# Plot metagene
#===================================================================================
#printf $lastname".solexa_result\t"$lastname1".solexa_result\t"$lastname > exp
#printf $lastname"_metagene_tss\n"$3".chrom.sizes\n" > exp2
#R --vanilla < metagene.r

#====================================================================================
#Run MACS
macs14 -t $lastname".bed"  -c $lastname1".bed" -g $4 --keep-dup=1 -n $lastname"_10-9" -p 1e-9
macs14 -t $lastname".bed"  -c $lastname1".bed" -g $4 --keep-dup=1 -n $lastname"_10-7" -p 1e-7

#====================================================================================
# Supplementary table1 and table2
#====================================================================================
if [[ "$3" = "hg18" ]]
then
	ref_file=`echo refflat_2006.txt`
elif [[ "$3" = "hg19" ]]
then
	ref_file=`echo refFlat_Hg19_Jan2014_fixeddates.txt`
elif [[  "$3" = "mm9" ]]
then
	ref_file=`echo refflat_mm9.txt`
fi
perl tss_sorted.pl $ref_file $lastname"_10-7"_summits.bed > $lastname"_10-7".supp2.txt
perl tss_sorted.pl $ref_file $lastname"_10-9"_summits.bed > $lastname"_10-9".supp2.txt

perl tss_closest.pl $ref_file $lastname"_10-7"_summits.bed > $lastname"_10-7".supp1.txt
perl tss_closest.pl $ref_file $lastname"_10-9"_summits.bed > $lastname"_10-9".supp1.txt


