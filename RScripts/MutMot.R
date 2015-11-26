#!/usr/bin/Rscript

#devtools::install_github("hadley/readr")
# usage 1. samtools view Myb.bam > Myb.sam
# usage 2. nohup./MutMot.R Myb.sam ACCGTT hg19 8 &
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(dplyr)
library(readr)

args <-commandArgs(trailingOnly=TRUE)
SAMFILE <-args[1]
MOTIF <- args[2]
GENOME <- args[3]
GENOME <- BSgenome::getBSgenome(GENOME)
DEPTH <-args[4]
DEPTHFRAC <- args[5]
REVMOTIF <- as.character(reverseComplement(DNAString(MOTIF)))
PATTERN <- paste0(MOTIF, "|", REVMOTIF)
LENMOTIF<-nchar(MOTIF)

SAM<-read_tsv(SAMFILE, col_names=F)
SAM<- SAM %>%
  mutate(ReadLength = nchar(X10))

gr<- GRanges(seqnames=Rle(SAM$X3),
             ranges=IRanges(SAM$X4, end=SAM$X4 + SAM$ReadLength))

gr.cov <- coverage(gr)
rm(gr)
SAM<- SAM %>% filter(grepl(PATTERN, X10))

SAM<- SAM %>% 
  mutate(RefSeq = as.character(BSgenome::getSeq(GENOME, X3, start= X4, end= X4  + ReadLength))) %>%
  filter(!grepl(PATTERN, RefSeq)) %>%
  mutate(pos = regexpr(PATTERN, X10), MotifStart=X4+pos) %>%
  group_by(X3, MotifStart) %>%
  mutate(MotifDepth=n()) %>%
  ungroup() %>%
  filter(MotifDepth > DEPTH)

 
SAMSUM <- group_by(SAM, X3, MotifStart) %>%
  summarise(MotifDepth=n()) %>%
  mutate(MotifEnd = MotifStart + LENMOTIF)


SAMSUM$AllDepth <- colMeans(apply(SAMSUM, 1, function(MOTIF) as.vector(gr.cov[[MOTIF['X3']]])[MOTIF['MotifStart']:MOTIF['MotifEnd']]))

SAMSUM <- mutate(SAMSUM, MotifFraction=MotifDepth/AllDepth) %>%
  ungroup() %>%
  arrange(desc(MotifDepth)) %>%
  select(X3, MotifStart, MotifEnd, everything()) %>%
  mutate(SampleID=gsub(".sam", "", SAMFILE)) %>% 
  filter(MotifFraction > DEPTHFRAC)

colnames(SAMSUM)[1]<-"Chr"
write_tsv(SAMSUM, paste0(MOTIF, ".", gsub(".sam", "", SAMFILE),  ".tsv"))

