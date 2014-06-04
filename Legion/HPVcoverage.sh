
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
BT2GENOME="/home/rmgzshd/Scratch/Genome/HPV_BT2/HPV"
SAMBAMBA="/home/rmgzshd/Tools/sambamba/sambamba"
GENOMESIZES="/home/rmgzshd/Scratch/Genome/HPV_BT2/HPV.sizes"
LOGFILE="/home/rmgzshd/Scratch/HNSC/Results/coverage.log"
RESULTS="/home/rmgzshd/Scratch/HNSC/Results"
WD="/home/rmgzshd/Scratch/HNSC"
export BT2GENOME; export SAMBAMBA; export RESULTS; export LOGFILE; export WD;

# 8. Your work *must* be done in $TMPDIR 
#cd $TMPDIR

FOUND=$(find -maxdepth 2 -name *.sorted.bam)
for BAM in $FOUND
do
	
	PREFIX=$(echo ${BAM} | sed 's/.sorted.bam//')
	printf "\n\n beginning %s sample.\n " $PREFIX
	
	$SAMBAMBA view --format=bam -F  "paired and (unmapped or mate_is_unmapped)" $BAM 				|
	tee bedtools bamtofastq -i stdin -fq R1.fq -fq2 R2.fq2										|
	bedtools bamtobed -bedpe -i stdin 	| sort -S8G -k 7											|
	awk '{if ($2 != -1 || $5 != -1) print $1,$2,$3,$4,$5,$6, $7, $8, $9, $10}' > $PREFIX.anchor.bedpe
			
	# Use bowtie2 to find matches amongst the unmapped genome in HPV
	# convert to bam
	# sort it 	
	bowtie2 -p 24 -x $BT2GENOME -1 R1.fq -2 R2.fq 				|
	$SAMBAMBA view -S --format=bam /dev/stdin > tempsortbam	 
	samtools sort -@ 12 -m 32G tempsortbam $PREFIX.hpv
	
	# NB sambabmba sort whilst fast doesn't work on stream???
	# So it fucks up the pipe flow here (see above)
	bedtools bamtobed -bedpe -i $PREFIX.hpv.bam | sort -S8G -k 7 			|
	awk '{if ($2 != -1 || $5 != -1) print $1,$2,$3,$4,$5, $6, $7, $8, $9, $10}' > $PREFIX.hpv.bedpe

    # calculate coverage NB NEED TO WRITE/APPEND THIS TO A LOG FILE
    printf "\n HPV coverage :"
    HPVCOV=$(bedtools genomecov -bg -ibam $PREFIX.hpv.bam  -g $GENOMESIZES  |
   	awk '{sum+=  $4} END { print sum/NR}')
	printf "HPV\t$PREFIX\t$HPVCOV\n" >> $LOGFILE

	# merge the single mapped HPV and ANCHOR samples to look for bridges
	join -1 7 -2 7 $PREFIX.hpv.bedpe $PREFIX.anchor.bedpe | 
	sort -k1,14 -k2,15 >  $PREFIX.integration.bed

	rm $PREFIX.hpv.bedpe 
	rm $PREFIX.anchor.bedpe
	mv $PREFIX.hpv.bam $RESULTS
	mv $PREFIX.integration.bed $RESULTS
	printf "\n Next file?: "
done
printf "\n No, all done! \n"


rm R1.fq
rm R2.fq
rm tempsortbam
