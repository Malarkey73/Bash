#!/bin/bash -l
#$ -S /bin/bash
# 2. Request ten hours of wallclock time (format hours:minutes:seconds).
#$ -l h_rt=10:00:0
# 3. Request 1 gigabyte of RAM 
#$ -l mem=32G
# 4. Request 15 gigabyte of TMPDIR space (default is 10 GB)
#$ -l tmpfs=32G
# 5. Set the name of the job.
#$ -N STAR
# 6. Find <your_project_id> by running the command "groups"
#$ -P TCGAHPVHNSC
# 7. Replace "<your_UCL_id>" with your UCL user ID :)
#$ -wd /home/rmgzshd/Scratch/Pablo/FASTQ
# 8. request 16 threads
#$ -l thr=16

set -o nounset
set -o errexit
set -o pipefail

ENSEMBL75="/home/rmgzshd/Scratch/Genome/ENS_hs_75_STARindex"
SEQTK="/home/rmgzshd/Tools/seqtk/seqtk"
module load star/2.3.0e
module load samtools/0.1.19
export ENSEMBL75; export FQDATA; export SEQTK;

for fq1 in *_R1_001.fastq.gz
do
    prefix=$(echo ${fq1} | sed 's/_R1_001.fastq.gz//')
    fq2=${prefix}_R2_001.fastq.gz

    STAR --genomeDir $ENSEMBL75 \
	--readFilesIn $fq1 $fq2 \
	--runThreadN 16 \
	--readFilesCommand $SEQTK trimfq -b 5 -e 10 \
	--genomeLoad NoSharedMemory \
	--outFileNamePrefix $prefix \
	--outSAMstrandField intronMotif \
	--outStd SAM \
	--outSAMattributes Standard | samtools view -buS - | samtools sort -m16G - $prefix

done


for bam in *.bam
do
	prefix=$(echo ${bam} | sed 's/.bam//')
	samtools index $bam $prefix.bam.bai
done
