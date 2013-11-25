#!/bin/bash

# data shortcuts
FQFOLDER="/mnt/store1/rawdata/FASTQ/dogliver"
BAMFOLDER="/mnt/store1/rawdata/FASTQ/dogliver/BAM"
# reference shortcuts
GENOME="/mnt/store1/Canis_familiaris/UCSC/canFam3/Sequence/Bowtie2Index/genome"
GENES="/mnt/store1/Canis_familiaris/UCSC/canFam3/Annotation/Genes/genes.gtf"
# tool shortcuts - NB bowtie needs to be in PATH because it isn't specified in script
SEQTK="/home/rmgzshd/seqtk/seqtk"
TOPHAT2="/home/rmgzshd/tophat2/tophat2"

# samtools and bowtie2 are already in my path
export FQFOLDER; export BAMFOLDER; export GENOME; export GENES; export SEQTK; export TOPHAT2;

# This is a basic script to trim, align, convert bam to sam and then sort unpaired reads. To adapt it to you should:

# 1. Put it in the top level directory with your fastq files (same as FQFOLDER above) . Whilst not strictly necesary...
# this is good scientific practice to keep a log WITH the data and all derived data in folders below. All other schemes
# eventually end in confusion and angst.

# 2. Change the shortcuts to reflect your own system paths to the data, reference genomes (mm9, mm10 etc) and tools

# 3. Check the name of your files and change the for loop and prefix string (line 27 & 30) accordingly.

#loop through fastq files
for fq in $FQFOLDER/*.fastq
do
        # get the sample name prefix and make a folder for output of each fastq
        prefix=$(echo ${fq} | sed 's/.fastq//')
        echo $prefix
        mkdir -p $prefix
        echo "making " $prefix

        # create a fifo for trimmed reads
        if [[ ! -p $FQFOLDER/R1.fifo ]]; then
                mkfifo $FQFOLDER/R1.fifo
        fi

        # trim/pipe into tophat2
        $SEQTK trimfq $fq > R1.fifo & \
        # -p 24 = 24 threads, again piped directly into samtools to convert to bam sorted with no intermediate files
        tophat2 -p 24 -G $GENES -o $prefix $GENOME $FQFOLDER/R1.fifo
done

# tidy up; deleting the temporary fifo
rm *.fifo