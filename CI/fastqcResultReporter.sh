#!/bin/bash

# given a folder full of fastqc results from multiple samples this script will
# make a simpel summary table of all those results in a single file
# you might use it a s part of a QC process with the similar flagstatResultParser.sh
# which is for BAM file stats.

FQCFILES=`find . -print | grep -i 'fastqc_data.txt'`
FIRST=true

for FQC in $FQCFILES
do
	# Samplename
	grep -n "Filename" $FQC | cut -f 2 | tr "\n" "\t" >> fastqc.data
	#Total Sequences
	grep -n "Total Sequences" $FQC | cut -f 2 | tr "\n" "\t" >> fastqc.data
	# GC percent
	grep -n -m 1 "%GC" $FQC | cut -f 2 | tr "\n" "\t" >> fastqc.data
	printf "\n"  >>fastqc.data

	if [[ $FIRST = true ]]; then
		sed -n '/>>Per base sequence quality\t/{:a;n;/>>END_MODULE/b;p;ba}' $FQC | cut -f 1,2 > pb.seq.quality
		
		sed -n '/>>Per sequence quality scores\t/{:a;n;/>>END_MODULE/b;p;ba}' $FQC | cut -f 1,2 > ps.quality.scores
		
		sed -n '/>>Per base GC content\t/{:a;n;/>>END_MODULE/b;p;ba}' $FQC | cut -f 1,2 > pb.seq.GC
		
		sed -n '/>>Per base N content\t/{:a;n;/>>END_MODULE/b;p;ba}' $FQC | cut -f 1,2 > pb.seq.N

		sed -n '/>>Sequence Duplication Levels\t/{:a;n;/>>END_MODULE/b;p;ba}' $FQC | cut -f 1,2 | tail -n +2 > dup.seq.lev

		sed -n '/>>Kmer Content\t/{:a;n;/>>END_MODULE/b;p;ba}' $FQC | cut -f 1,2,3 > kmer.lev

		sed -n '/>>Per base sequence content\t/{:a;n;/>>END_MODULE/b;p;ba}' $FQC | cut -f 1,2 > tempfile
		sed -n '/>>Per base sequence content\t/{:a;n;/>>END_MODULE/b;p;ba}' $FQC | cut -f 1,3 | tail -n +2 >> tempfile
		sed -n '/>>Per base sequence content\t/{:a;n;/>>END_MODULE/b;p;ba}' $FQC | cut -f 1,4 | tail -n +2 >> tempfile
		sed -n '/>>Per base sequence content\t/{:a;n;/>>END_MODULE/b;p;ba}' $FQC | cut -f 1,5 | tail -n +2 >> tempfile
		cat tempfile > pb.sequence.scores

		FIRST=false
	else 
		paste pb.meanseq.quality <(sed -n '/>>Per base sequence quality\t/{:a;n;/>>END_MODULE/b;p;ba}' $FQC | cut -f 2) > tempfile
		cat tempfile > pb.seq.quality
		
		paste ps.quality.scores <(sed -n '/>>Per sequence quality scores\t/{:a;n;/>>END_MODULE/b;p;ba}' $FQC | cut -f 2) > tempfile
		cat tempfile > ps.quality.scores
		
		paste pb.seq.GC <(sed -n '/>>Per base GC content\t/{:a;n;/>>END_MODULE/b;p;ba}' $FQC | cut -f 2) > tempfile
		cat tempfile > pb.seq.GC
		
		paste pb.seq.N <(sed -n '/>>Per base N content\t/{:a;n;/>>END_MODULE/b;p;ba}' $FQC | cut -f 2) > tempfile
		cat tempfile > pb.seq.N

		paste dup.seq.lev <(sed -n '/>>Sequence Duplication Levels\t/{:a;n;/>>END_MODULE/b;p;ba}' $FQC | cut -f 2 | tail -n +2) > tempfile
		cat tempfile > dup.seq.lev

		paste kmer.lev <(sed -n '/>>Kmer Content\t/{:a;n;/>>END_MODULE/b;p;ba}' $FQC | cut -f 1,2,3) > tempfile
		cat tempfile > kmer.lev

		sed -n '/>>Per base sequence content\t/{:a;n;/>>END_MODULE/b;p;ba}' $FQC | cut -f 2 > tempfile
		sed -n '/>>Per base sequence content\t/{:a;n;/>>END_MODULE/b;p;ba}' $FQC | cut -f 3 | tail -n +2 >> tempfile
		sed -n '/>>Per base sequence content\t/{:a;n;/>>END_MODULE/b;p;ba}' $FQC | cut -f 4 | tail -n +2 >> tempfile
		sed -n '/>>Per base sequence content\t/{:a;n;/>>END_MODULE/b;p;ba}' $FQC | cut -f 5 | tail -n +2 >> tempfile
		paste pb.sequence.scores tempfile > tempfile2
		cat tempfile2 > pb.sequence.scores

	fi
done
