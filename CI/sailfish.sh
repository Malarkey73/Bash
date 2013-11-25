#!/bin/bash

# a sailfish gene index first needs built, e.g. something like:
#nohup sailfish index -t /mnt/store1/Mus_musculus/gencode.vM1.pc_transcripts.fa -o /mnt/store1/Mus_musculus/sailfishIndex -k 20 &

#

#genome shortcuts
SFINDEX="/mnt/store1/Mus_musculus/sailfishIndex"

# data shortcuts
FQFOLDER="/mnt/store1/rawdata/FASTQ/mouseliver"

# software
SAILFISH="/home/rmgzshd/Sailfish/bin/sailfish"
SEQTK="/home/rmgzshd/seqtk/seqtk"


export SFINDEX; export FQFOLDER; export SAILFISH

for fq in $FQFOLDER/*.fastq
do
	# get the sample name prefix
        prefix=$(echo ${fq} | sed 's/.fastq//')
    $SEQTK trimfq $fq > $prefix.trim    
done
        
$SAILFISH quant -p 24 -i $SFINDEX -r $FQFOLDER/*.trim -o $FQFOLDER




