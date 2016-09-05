### Precepts
Before writing pipeline scripts I begin with a set of rules about 1) file organisation, 2) modular outputs, and 3) naming conventions. By sticking to these rules I can make the scripts simpler and similar.

## 1) File Organisation
The starting point for most genomic data is a colection of raw fastq files. For each project I begin with a Project folder that contains a FASTQ folder and I copy or download the QC, Trimming or, Alignment script to that folder. If I run a FastQC script this produces results (zip files) in a FASTQC folder within that. If alignment then the results go to a new BAM file folder within the FASTQ folder. If I were to then do a peak calling I would copy/load the script to the BAM file folder and results would be created in a BED file folder within that. So each folder - usually named after the result type within it also contains a copy of the script that produces results in the folders within it. I keep these tiny files as up to date permenant records of the analyses so far. Folders are uppercase, Files are lower or often Camel.Case - this makes scripts easier to understand.

## 2) Modular Outputs
I try to write scripts according to a unix(ish) philosophy - they do one thing well. This means I dont write scripts that convert FASTQ to BED results I write FASTQ to BAM, and BAM to BED. And each script begins with a defined file format and outputs a defined file format (BAM, BED, VCF, TSV, FASTQ). In one sense this is less automated - you probably have to check your intermediate results (no bad thing maybe). But it also makes pipelines composable you can mix various existing scripts or  - since inputs and outputs are defined standards, drop in your own script in whatever language you like.

## 3) Naming Conventions
I use naming conventions that carry information about the sample and file types downstream through the pipeline. This makes the code and usage of scripts much easier as the scripts will automatically recognise certain types of file and process them accordingly (e.g. paired, compressed, see below). Plus consistency will also help you to understand a project at a glance should you return to it after several months. So for instance examine:

ProjSH_S01.N.R2_001.fastq.gz

Each important part of the name is separated by ".", as Bash works well with these separators. (e.g. PREFIX="ProjSH_S01" ; mv $PREFIX.* $HOME)
This is sample 1 (S01) from the ProjSH project and folder. This prefix will follow the pipeline downstream from beginning to end. You will find it very useful if this prefix is alphanumerically in the same order as the samples - as software such as R or Unix will order them that way (there is nothing more annoying than a graph/plot with sample 1,10,11, then 2,3,4 etc). The "N" here denotes it is a normal sample, whereas a tumour sample would be "T". Then "R2_001" is an Illumina convention for the second of paired end files (scripts will throw errors if there is not a matching "R1_001" with "R2_001"). Next is the file type "fastq" and then perhaps a "gz" suffix denoting compression.

### Script Design
I write bash scripts - but you should be able to drop in your own to work with mine in any language that the CS server will handle. The most important script I write is called "qsubber.sh". This takes as input the type of job you want to run and analyses the filenames above to extract the information within them (paired, compressed, tumour/normal, etc) and then for each file or pair of files (e.g. R1_001 an dmatching R2_001) it calls a new qsub job - called with argument parameters pasing that information on. So all the difficult regular expression code is in that one wrapper script, and all jobs are parallelised. The other scripts are simpler and contain no loops and few logical switches. So for instance to do an alignment on the above file "ProjSH_S01.N.R2_001.fastq.gz" you might call qsubber like so.

qsubber.sh bwa.qsub

This would then produce a set of calls including for that file:

qsub bwa.qsub ProjSH_S01 1 1

where we hve the qsub command, the qsub script, the file prefix, and then two logical flags indicate paired and compressed.








