#!/bin/bash

SFINDEX="/mnt/store1/Homo_sapiens/Gencode19/sailfish_gene"
OUTDIR="/mnt/store1/sailfish_results"

export SFINDEX;  export OUTDIR;

for fq in *_R1_001.trimmed.fastq.gz
do
  	# get the sample name prefix
        prefix=$(echo ${fq} | sed 's/_R1_001.trimmed.fastq.gz//')

        # create a fifo for both R1 and R2 pairs
        if [[ ! -p MATE1.fifo ]]; then
                mkfifo MATE1.fifo
        fi

                if [[ ! -p MATE2.fifo ]]; then
                mkfifo MATE2.fifo
        fi

	#decompress the mate1 and 2 into a fifo
        zcat -c $fq > MATE1.fifo & \
        zcat -c ${prefix}_R2_001.fastq.gz > MATE2.fifo &

        # and stream the fifo into sailfish
        sailfish quant -p 24 -i $SFINDEX -l "T=PE:O=><:S=U" \
        -1 MATE1.fifo -2 MATE2.fifo -o $OUTDIR/$prefix

done