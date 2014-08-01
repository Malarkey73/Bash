#!/bin/bash -l
#$ -S /bin/bash
# 2. Request ten hours of wallclock time (format hours:minutes:seconds).
#$ -l h_rt=10:00:0
# 3. Request 1 gigabyte of RAM 
#$ -l mem=32G
# 4. Request 15 gigabyte of TMPDIR space (default is 10 GB)
#$ -l tmpfs=32G
# 5. Set the name of the job.
#$ -N FAIDX
# 6. Find <your_project_id> by running the command "groups"
#$ -P TCGAHPVHNSC
# 7. Replace "<your_UCL_id>" with your UCL user ID :)
#$ -wd /home/rmgzshd/Scratch/CAKE/BAM
# 8. request 16 threads
#$ -l thr=16

set -o nounset
set -o errexit
set -o pipefail


module load samtools/0.1.19
samtools faidx /home/rmgzshd/Scratch/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa

