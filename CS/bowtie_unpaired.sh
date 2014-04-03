#!/bin/bash

# data shortcuts
FQFOLDER="~/testfq"
BAMFOLDER="~/testfq/bam"
# reference shortcuts
GENOME="/SAN/biomed/biomed13/cohesin-bio/TOOL/Cufflinks/GTF/Mus_musculus/Ensembl/NCBIM37/Sequence/Bowtie2Index/genome"
GENES="/SAN/biomed/biomed13/cohesin-bio/TOOL/Cufflinks/GTF/Mus_musculus/Ensembl/NCBIM37/Annotation/Genes/genes.gtf"
# tool shortcuts
SEQTK="/SAN/biomed/biomed13/cohesin-bio/TOOL/seqtk/seqtk"
BOWTIE2="/SAN/biomed/biomed13/cohesin-bio/TOOL/Bowtie2/bowtie2-2.1.0/bowtie2"
SAMTOOLS="/SAN/biomed/biomed13/cohesin-bio/TOOL/SAMtools/samtools-0.1.19/samtools"

export FQFOLDER; export BAMFOLDER; export GENOME; export GENES; export SEQTK; export BOWTIE2; export SAMTOOLS

# This is a basic script to trim, align, convert bam to sam and then sort paired reads. To adapt it to you should:
# 1. Put it in the top level directory with your fastq files (same as FQFOLDER above).
# 2. Change the shortcuts to reflect your own system paths to the data, reference genomes (mm9, mm10 etc) and tools
# 3. Check the name of your files and change the for loop and prefix string (line 24 & 27) accordingly.

for fq in $FQFOLDER/*fastq
do
        # get the sample name prefix
        prefix=$(echo ${fq} | sed 's/fastq//')
        # create a fifo 
        if [[ ! -p $FQFOLDER/R1.fifo ]]; then
                mkfifo $FQFOLDER/R1.fifo
        fi
        # decompress and trim fifo on the run and pipe into bowtie
        $SEQTK trimfq $fq > R1.fifo &
        # -p 24 = 24 threads, again piped directly into samtools to convert to bam sorted with no intermediate files
        $BOWTIE2 -p 10 -x $GENOME -U $FQFOLDER/R1.fifo | $SAMTOOLS view -bS - | $SAMTOOLS sort - $prefix".sorted"
done

# move the bam files to our bamfolder
mkdir -p $BAMFOLDER
mv $FQFOLDER/*.bam $BAMFOLDER
# tidy up; deleting the temporary fifo
rm *.fifo
