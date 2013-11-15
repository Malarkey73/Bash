#$ -S /bin/bash
#$ -l h_vmem=10G
#$ -l tmem=10G
#$ -l h_rt=1000:0:0
#$ -cwd
#$ -j y
#$ -V
  
#shortcuts to tools
FASTQC="/SAN/biomed/biomed13/cohesin-bio/TOOL/FastQC/FastQC/fastqc"
# user data - note sure these are good places
FQFOLDER="/SAN/biomed/biomed13/cohesin-bio/RNASEQ/result/20120423/FASTQ/"

export JAVA_HOME="/share/apps/jdk1.7.0_45"
export PATH=$JAVA_HOME/bin:$PATH

# make sure child processes see shortcuts too
export FASTQC; export FQFOLDER;

# c.a. 5 mins per file (it only analyses first 100-200k reads I think)
for fq in $FQFOLDER/*.fq.gz
do
# a report directory for each file ending .fq.gz every 
  	$FASTQC  $fq 
done
rm $FQFOLDER/*.zip






