
#!/bin/bash

SEQTK="/home/rmgzshd/seqtk/seqtk"
STAR="/home/rmgzshd/STAR/bin/Linux_x86_64/STAR"
GTDOWNLOAD="/home/rmgzshd/cghub/bin/gtdownload"
CGHUBKEY="/home/rmgzshd/cghub/cghub.key"


export SEQTK; export STAR; export GTDOWNLOAD; export CGHUBKEY;

# script args 
# e.g. CESC1.tsv
IDSTSV=`echo $1`
# usually column 23
IDSCOLUMN=`echo $2`

# reads the IDS of selected column from the file into IDSARRAY, skip header (NR>1)
IDSARRAY=$(awk < $IDSTSV -v x=$IDSCOLUMN 'NR >1{ print $x}')

# is this doing the loop for every line?
for ID in $IDSARRAY
#for ID in "${IDSARRAY[@]}"
do
	# download the tar.gz file from UCSC CGHUb
   	$GTDOWNLOAD --max-children 12 -c $CGHUBKEY $ID
   	wait
   	# remove the temporary download file
   	rm "$ID.gto"
   	# gunzip the archive to reveal 2 PE fastq files 
   	tar -xzvf $ID/*.tar.gz -C $ID
   	# move them to their folder
   	#mv *.fastq $ID
   	# capture the names and locationsof the PE fq files
   	fq1=($ID/*_1.fastq)
   	fq2=($ID/*_2.fastq)

   	# run STAR aligner CESC RNA against the HPV genomes
   	$STAR --genomeDir /mnt/store1/CESC_RNA/STAR_HPV_GENOME \
	--readFilesIn $fq1 $fq2 \
	--runThreadN 16 \
	--readFilesCommand $SEQTK trimfq -b 5 -e 5 \
	--genomeLoad NoSharedMemory \
	--outSAMstrandField intronMotif \
	--outSAMtype BAM Unsorted \
	--outFileNamePrefix $ID/$ID. \
	--outSAMattributes Standard

	#remove used raw data as I have no space
	rm "$ID/*.fastq"
	rm "$ID/*.tar.gz"
done


