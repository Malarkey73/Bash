#!/usr/bin/Rscript


# nohup ./indelfixer.R NB2.indel.Somatic.hc

args<-commandArgs(trailingOnly=TRUE)
file<-args[1]
indel <- read.delim(file, as.is=T)

indel2 <-indel


ins <- grep("+", indel2$var, fixed=T)
del <- grep("-", indel2$var, fixed=T)


indel2$ref[ins] <- "-"
indel2$var[ins] <- gsub("+", "", indel2$var[ins], fixed=T)

indel2$ref[del]<- paste0(indel2$ref[del], gsub("-", "", indel2$var[del], fixed=T))
indel2$var[del]<- "-"

indel2$end <-ifelse(indel2$ref=="-", indel2$position, indel2$position + (nchar(indel2$ref)-1))

indelout<- cbind(indel2[,1:2], indel2$end, indel2[,3:4])
file <- gsub("hc", "pass", file)
write.table(indelout, file, col.names=FALSE, row.names=FALSE, quote=FALSE)
