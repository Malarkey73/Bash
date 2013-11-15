#!/bin/bash

# data shortcuts
FQFOLDER="/mnt/store1/rawdata/fastq/Suzana"
BAMFOLDER="/mnt/store1/rawdata/fastq/Suzana/bam"
# reference shortcuts
GENOME="/mnt/store1/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome"
GENES="/mnt/store1/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf"
# tool shortcuts
SEQTK="/home/rmgzshd/seqtk/seqtk"
BOWTIE2="/home/rmgzshd/bowtie2-2.1.0/bowtie2"
SAMTOOLS="/home/rmgzshd/samtools/samtools"

export FQFOLDER; export BAMFOLDER; export GENOME; export GENES; export SEQTK; export BOWTIE2; export SAMTOOLS

# This is a basic script to trim, align, convert bam to sam and then sort paired reads. To adapt it to you should:

# 1. Put it in the top level directory with your fastq files (same as FQFOLDER above) . Whilst not strictly necesary...
# this is good scientific practice to keep a log WITH the data and all derived data in folders below. All other schemes
# eventually end in confusion and angst.

# 2. Change the shortcuts to reflect your own system paths to the data, reference genomes (mm9, mm10 etc) and tools

# 3. Check the name of your files and change the for loop and prefix string (line 24 & 27) accordingly.

#loop only through the first pair
for fq in $FQFOLDER/*R1_001.fq.gz
do
        # get the sample name prefix
        prefix=$(echo ${fq} | sed 's/R1_001.fq.gz//')
        # create a fifo for both R1 and R2 pairs
        if [[ ! -p $FQFOLDER/R1.fifo ]]; then
                mkfifo $FQFOLDER/R1.fifo
        fi

if [[ ! -p $FQFOLDER/R2.fifo ]]; then
                mkfifo $FQFOLDER/R2.fifo
        fi

        # decompress and trim each fifo pair on the run and pipe into bowtie
        $SEQTK trimfq $fq > R1.fifo & \
        $SEQTK trimfq ${prefix}R2_001.fq.gz > R2.fifo & \
        # -p 24 = 24 threads, again piped directly into samtools to convert to bam sorted with no intermediate files
        $BOWTIE2 -p 24 -x $GENOME -1 $FQFOLDER/R1.fifo -2 $FQFOLDER/R2.fifo | $SAMTOOLS view -bS - | $SAMTOOLS sort - $prefix"sorted"
done

# move the bam files to our bamfolder
mkdir -p $BAMFOLDER
mv $FQFOLDER/*.bam $BAMFOLDER
# tidy up; deleting the temporary fifo
rm *.fifo