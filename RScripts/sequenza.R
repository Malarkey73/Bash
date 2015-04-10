#!/usr/bin/RScript

args<-commandArgs(TRUE)

snp <- read.table(args[1], header = TRUE, sep = "\t")
cnv <- read.table(args[2], header = TRUE, sep = "\t")
seqz.data <- VarScan2seqz(varscan.somatic = snp, varscan.copynumber = cnv)
out=paste(args[1], ".seqz")
write.table(seqz.data, paste0(args[1], ".seqz"), col.names = TRUE, row.names = FALSE, sep = "\t")
