#!/usr/bin/env Rscript
# usage Rscript --vanilla ./MotifFinder2.R ACCGTT DU528_MYB.hg19_multianno.vcf DU528_MYB.bam
# usage Rscript --vanilla ./MotifFinder2.R ACCGTT PF382_MYB.hg19_multianno.vcf PF382_MYB.bam


# Currently used within scapel_pipeline.sh script.

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 3) {
  stop("Command line arguments must be supplied : Motif, VCF, and BAM - in that order", call.=FALSE)
}

library(Rsamtools)
library(VariantAnnotation)
library(Biostrings)
library(dplyr)
library(GenomicRanges)
library(readr)

# take the motif from args, and join with the reverse complement into a single pattern search string joined by |
Motif<- args[1]
RevCMotif <- as.character(reverseComplement(DNAString(Motif)))
DualPattern <- paste0(Motif, "|", RevCMotif)

# args for name of VCF, and BAM and optionally the ref genome
VCF.file <- args[2]
BAM.file <- args[3]
OUT.var.file <- gsub(".bam", ".varMotifs.txt", BAM.file)
OUT.som.file <- gsub(".bam", ".varSomatic.txt", BAM.file)

refgenome <- ifelse(length(args) == 4, args[4], "hg19")
#if (length(args) == 4) {
#  refgenome <- args[4]
#}else{
#  refgenome <-"hg19"
#}


# The VCF file has already identified potential variants and also contains dbSNP and 1000 genome info
# It also contains some depth of coverage information
vcf<- readVcf(VCF.file, genome = refgenome)

# the following function seems to be renamed in later versions of VariantAnnotation
# VCFpos <- rowRanges(vcf)
VCFpos   <- rowData(vcf)

# Only import Bam reads that overlap potential variants in the VCF file
param <- ScanBamParam(which = VCFpos, what= c("rname","strand", "pos", "seq"))
bam<-scanBam(file=BAM.file, param=param)

# I put those reads into an easier to parse DataFrame format
lst <- lapply(names(bam[[1]]), function(elt) {
       do.call(c, unname(lapply(bam, "[[", elt)))
  })
names(lst) <- names(bam[[1]])
df <- do.call("DataFrame", lst)


# I work out which of the imported reads have a Motif site or its reverse complement (DualPattern)
# And I find the exact spot where it is (Motifpos).
# I then count the occurences of each site and filter to only those with more than 10 reads evidence
df2<- df %>%
  as.data.frame() %>%
  filter(grepl(DualPattern , seq)) %>%
  mutate(patpos = regexpr(DualPattern,seq), Motifpos= (patpos -1) + pos) %>%
  group_by(rname,Motifpos) %>%
  mutate(n=n()) %>%
  filter(n > 10)

# I take those positions and create a genomic range object
gr<- GRanges(seqnames=Rle(df2$rname),
             ranges=IRanges(start = df2$Motifpos-2, end = df2$Motifpos+ 8))

# From the reads in motif ranges I identify VCFS near the motif sites
VCFmotifs<- VCFpos[overlapsAny(VCFpos, gr),]
# Whereas from the Annovar dbsnp info I identify somatic mutations
VCFsomatic <- VCFpos[unlist(info(vcf)$snp138)==".",]

# You have to change ALT to character otherwise it doesn't export correctly!!!
VCFmotifs$ALT<-as.character(unlist(VCFmotifs$ALT))
VCFsomatic$ALT<-as.character(unlist(VCFsomatic$ALT))


mnmotifs<- match(names(VCFmotifs), names(VCFpos))
mnsomatic<- match(names(VCFsomatic), names(VCFpos))

# I combine the positional information from VCF (REF, ALT, etc) with the extra dbSNP, 1000g info
# Plus I add acolumn for the sampleID in case I want to merge all outputs into a single file later.
# [,-14] is a hacky hack to remove a column I hate
outMotifs <- cbind(DataFrame(VCFmotifs), info(vcf[mnmotifs,])[,-14])
outMotifs$sampleID <- gsub(".bam", "", BAM.file)
write_tsv(as.data.frame(outMotifs), OUT.var.file)

outSomatic <- cbind(DataFrame(VCFsomatic), info(vcf[mnsomatic,])[,-14])
outSomatic$sampleID <- gsub(".bam", "", BAM.file)
write_tsv(as.data.frame(outSomatic), OUT.som.file)



