
#!/bin/bash -1
#$ -S /bin/bash 

# 2. Request ten hours of wallclock time (format hours:minutes:seconds).
#$ -l h_rt=10:00:0

# 3. Request 1 gigabyte of RAM 
#$ -l mem=8G

# 4. Request 15 gigabyte of TMPDIR space (default is 10 GB)
#$ -l tmpfs=32G

# 5. Set the name of the job.
#$ -N HPVcoverage

# 6. Find <your_project_id> by running the command "groups"
#$ -P TCGAHPVHNSC

# 7. Replace "<your_UCL_id>" with your UCL user ID :)
#$ -wd /home/rmgzshd/Scratch/HNSC

set -o nounset
set -o errexit
set -o pipefail

#module unload compilers/intel/11.1/072
#module load compilers/gnu/4.6.3
#module load bowtie2
#module load samtools
#module load bedtools
