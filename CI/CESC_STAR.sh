
#!/bin/bash

SEQTK="/home/rmgzshd/seqtk/seqtk"
STAR="/home/rmgzshd/STAR/bin/Linux_x86_64/STAR"
GTDOWNLOAD="/home/rmgzshd/cghub/bin/gtdownload"
CGHUBKEY="home/rmgzshd/cghub/cghub.key"

fq1="140513_UNC15-SN850_0365_BC4BYLACXX_GTGAAA_L004_1.fastq"
fq2="140513_UNC15-SN850_0365_BC4BYLACXX_GTGAAA_L004_2.fastq"

export SEQTK; export STAR;

# e.g. CESC1.tsv
IDSTSV=`echo $1`
# usually column 22
IDSCOLUMN=`echo $2`

# reads the IDS of selected column from the file into IDSARRAY, skip header (NR>1)
IDSARRAY=$(awk < $IDSTSV -v x=$IDSCOLUMN 'NR >1{ print $x}')

for ID in "${IDSARRAY[@]}"
do
   $GTDOWNLOAD --max-children 12 -c $CGHUBKEY $ID


done


$STAR --genomeDir /mnt/store1/CESC_RNA/STAR_HPV_GENOME \
	--readFilesIn $fq1 $fq2 \
	--runThreadN 16 \
	--readFilesCommand $SEQTK trimfq -b 5 -e 10 \
	--genomeLoad NoSharedMemory \
	--outSAMstrandField intronMotif \
	--outSAMtype BAM Unsorted \
	--outFileNamePrefix TCGA-CESC-U2. \
	--outSAMattributes Standard &
