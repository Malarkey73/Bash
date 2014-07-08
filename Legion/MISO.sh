#!/bin/bash -l
#$ -S /bin/bash
# 2. Request ten hours of wallclock time (format hours:minutes:seconds).
#$ -l h_rt=5:00:0
# 3. Request 1 gigabyte of RAM 
#$ -l mem=32G
# 4. Request 15 gigabyte of TMPDIR space (default is 10 GB)
#$ -l tmpfs=32G
# 5. Set the name of the job.
#$ -N MISO
# 6. Find <your_project_id> by running the command "groups"
#$ -P TCGAHPVHNSC
# 7. Replace "<your_UCL_id>" with your UCL user ID :)
#$ -wd /home/rmgzshd/Scratch/Pablo/MISO
# 8. request 16 threads
#$ -l thr=16

set -o nounset
set -o errexit
set -o pipefail

BAMDATA="/home/rmgzshd/Scratch/Pablo/BAM"

module unload compilers/intel/11.1/072
module load compilers/gnu/4.6.3
module load python/enthought/7.3-2_2013-10-04
module load samtools/0.1.19
module load bedtools/2.17.0 
module load miso/0.5.2

export BAMDATA;

# paired end taken from header of "insert-dist/SC1_ATCACG_L006.bam.insert_len"
# after running MISOutils.sh script
miso --run ./indexed $BAMDATA/SC1_ATCACG_L006.bam \
--output-dir 	output 													\
--prefilter																\
--settings-filename /home/rmgzshd/Scratch/Pablo/MISO/settings/miso_settings.txt	\
--read-len 		100 													\
--paired-end 	174	53													\
--SGEarray																\
--job-name		MISO 													\
--use-cluster