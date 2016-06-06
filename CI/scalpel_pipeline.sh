#!/bin/bash

GENOME="/mnt/store1/GATK_BUNDLE_2.8_hg19/Homo_sapiens_assembly19.fasta"
SCALPEL="/home/rmgzshd/scalpel-0.4.1/scalpel"
TABLE_ANNOVAR="/home/rmgzshd/annovar/table_annovar.pl"
#annotations
DB="/home/rmgzshd/annovar/humandb"
#inputs
export GENOME; export SCALPEL; export TABLE_

# Prerequisites here are an exisiting BAM (using bwa with default settings) and BED file defining ChIP peaks.

for BAM in *.bam
do

  
	PREFIX=$(echo ${BAM} | sed 's/.bam//')
	BED="$PREFIX"_peaks.narrowPeak
	# this is the microassembly step within ChIP peak regions
	$SCALPEL --single --bam $BAM --bed $BED --numprocs 20 --ref $GENOME --dir ./"$PREFIX"
	
	# this exports VCF results in these regions
	$SCALPEL --export --db $PREFIX/variants.db --type all --bed $BED --covratio 0.1 --mincov 8 --ref $GENOME --format vcf > $PREFIX/$PREFIX.vcf

  	# this adds varaint annotations to the VCF
	$TABLE_ANNOVAR $PREFIX/$PREFIX.vcf $DB \
	-buildver hg19 \
	-out $PREFIX \
	-protocol refGene,1000g2015aug_all,snp138,cosmic70 \
	-operation g,f,f,f \
	-nastring . \
	-vcfinput

	# but Annovar produces alot of un-neccessary files we can remove to save space
	rm $PREFIX/*multianno.txt
	rm $PREFIX/*_snp138_*
	rm $PREFIX/*refGene*
	rm $PREFIX/*.log
	rm $PREFIX/*_cosmic70_*
	rm $PREFIX/*hg19_ALL.*
	rm $PREFIX/*.avinput
  
  # This compares the VCF and the BAM to find putative alterations that create novel Myb sites.
  # Then it exports all the above info in a TSV table format.
  Rscript --vanilla ./MotifFinder2.R ACCGTT DU528_MYB.hg19_multianno.vcf DU528_MYB.bam

done
