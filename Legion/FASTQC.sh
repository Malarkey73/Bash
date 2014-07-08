#!/bin/bash -l
#$ -S /bin/bash
# 2. Request ten hours of wallclock time (format hours:minutes:seconds).
#$ -l h_rt=3:00:0
# 3. Request 1 gigabyte of RAM 
#$ -l mem=8G
# 4. Request 15 gigabyte of TMPDIR space (default is 10 GB)
#$ -l tmpfs=4G
# 5. Set the name of the job.
#$ -N FASTQC
# 6. Find <your_project_id> by running the command "groups"
#$ -P TCGAHPVHNSC
# 7. Replace "<your_UCL_id>" with your UCL user ID :)
#$ -wd /home/rmgzshd/Scratch/Pablo/FASTQ
# 8. request 16 threads
#$ -l thr=8

set -o nounset
set -o errexit
set -o pipefail

FQDATA="/home/rmgzshd/Scratch/Pablo/FASTQ"
module load java/1.6.0_32
module load fastqc/0.10.1
export FQDATA


for fq in *.fastq.gz
do
	fastqc $fq
done

