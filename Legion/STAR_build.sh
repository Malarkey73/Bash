#!/bin/bash -l
#$ -S /bin/bash 
# 2. Request ten hours of wallclock time (format hours:minutes:seconds).
#$ -l h_rt=5:00:0
# 3. Request 1 gigabyte of RAM 
#$ -l mem=20G
# 4. Request 15 gigabyte of TMPDIR space (default is 10 GB)
#$ -l tmpfs=32G
# 5. Set the name of the job.
#$ -N STARbuild
# 6. Find <your_project_id> by running the command "groups"
#$ -P TCGAHPVHNSC
#7. Select 12 threads (the most possible on Legion).
#$ -l thr=12
# 8. Replace "<your_UCL_id>" with your UCL user ID :)
#$ -wd /home/rmgzshd/Scratch/Homo_sapiens/UCSC/hg19

set -o nounset
set -o errexit
set -o pipefail

module load star/2.3.0e
# Ensembl style reference
STAR 	--runMode genomeGenerate \
		--runThreadN 12 \
		--genomeDir ./Sequence/WholeGenomeFasta \
		--genomeFastaFiles ./Sequence/WholeGenomeFasta/genome.fa \
		--sjdbGTFfile ./Annotation/Genes/genes.gtf
		--sjdbOverhang 100
