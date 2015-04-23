#!/usr/bin/Rscript

# usage 
# here the assumption is existing matched NB7.snp and NB7.copynumber files
# nohup ./sequenza.R NB7

args<-commandArgs(trailingOnly=TRUE)
prefix<-args[1]
#read and convert Varscan output
snp <- read.table(paste0(prefix,".snp"), header = TRUE, sep = "\t")
cnv <- read.table(paste0(prefix,".copynumber"), header = TRUE, sep = "\t")
seqz.data <- sequenza::VarScan2seqz(varscan.somatic = snp, varscan.copynumber = cnv)


# gc normalisation
gc.stats <- gc.norm(x = seqz.data$depth.ratio, gc = seqz.data$GC.percent)
gc.vect  <- setNames(gc.stats$raw.mean, gc.stats$gc.values)
seqz.data$adjusted.ratio <- seqz.data$depth.ratio/gc.vect[as.character(seqz.data$GC.percent)]

# save data
out <- paste0(prefix, ".seq.gz")
write.table(seqz.data, gzfile(out), col.names = TRUE, row.names = FALSE, sep = "\t")

extract <- sequenza.extract(out)
CP   <- sequenza.fit(extract)
sequenza.results(extract, CP, sample.id = prefix)
