#!/bin/bash

RESULTDIRS=`find . -maxdepth 1 -type d`

for RESULTDIR in $RESULTDIRS
do
	prefix=$(echo $RESULTDIR | sed 's/^..//')
	bam=$RESULTDIR/*.bam
	samtools view $bam > $prefix.sam
	sed -i "s/$/\t$prefix/" $prefix.sam;
done

ls *.sam > samlist
cat samlist | xargs cat >> HPV.sam
rm *.sam
rm samlist
