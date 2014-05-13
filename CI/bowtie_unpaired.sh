#!/bin/bash
# data shortcuts
FQFOLDER="/mnt/store1/memeCTCF"
BAMFOLDER="/mnt/store1/memeCTCF/BAM"
# reference shortcuts
#GENOME="/mnt/store1/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome"
GENOME="/mnt/store1/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome"
# tool shortcuts
SEQTK="/home/rmgzshd/seqtk/seqtk"
BOWTIE2="/home/rmgzshd/bowtie2-2.1.0/bowtie2"
SAMTOOLS="/home/rmgzshd/samtools/samtools"

export FQFOLDER; export BAMFOLDER; export GENOME; export SEQTK; export BOWTIE2; export SAMTOOLS

# This is a basic script to trim, align, convert bam to sam and then sort reads and remove duplicates. To adapt it to you should:
# 1. Put it in the top level directory with your fastq files (same as FQFOLDER above) . Whilst not strictly necesary...
# this is good scientific practice to keep a log WITH the data and all derived data in folders below. All other schemes
# eventually end in confusion and angst.
# 2. Change the shortcuts to reflect your own system paths to the data, reference genomes (mm9, mm10 etc) and tools
# 3. Check the name of your files and change the for loop and prefix string (line 26 & 29) accordingly.

#loop through all fastq
for fq in $FQFOLDER/*fastq.gz
do
        # get the sample name prefix
        prefix=$(echo ${fq} | sed 's/.fastq.gz//')
        # create a fifo 
        if [[ ! -p $FQFOLDER/R1.fifo ]]; then
                mkfifo $FQFOLDER/R1.fifo
        fi
        # decompress and trim fifo on the run and pipe into bowtie
        $SEQTK trimfq $fq > R1.fifo & \
        # -p 24 = 24 threads, again piped directly into samtools to convert to bam sorted with no intermediate files
        $BOWTIE2 -p 24 -N 2 -x $GENOME -U $FQFOLDER/R1.fifo | $SAMTOOLS view -Sbu - | 
        $SAMTOOLS sort -m 8000000000 - $prefix".bam"  
done

# move the bam files to our bamfolder
mkdir -p $BAMFOLDER
mv $FQFOLDER/*.bam $BAMFOLDER
# tidy up; deleting the temporary fifo
rm *.fifo