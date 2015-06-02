#!/usr/bin/Rscript

# nohup ./vscantosnnov.R NB2 5 0.05

library(dplyr)
args<-commandArgs(trailingOnly=TRUE)

prefix<-args[1]
tvf<- args[2]
pval=args[3]

indel=paste0(prefix, ".indel")
snp=paste0(prefix, ".snp")


indel.tab <- dplyr::as.tbl(read.delim(indel, as.is=T, header=T))
snp.tab <- dplyr::as.tbl(read.delim(snp, as.is=T, header=T))

snp.tab<- snp.tab %>%
  dplyr::mutate(tumor_var_freq = as.numeric(gsub("%", "", tumor_var_freq)), position.end = position) %>%
  dplyr::filter(tumor_var_freq > tvf, somatic_status =="Somatic", somatic_p_value < pval) 
  
indel.tab<- indel.tab %>%
  dplyr::mutate(tumor_var_freq = as.numeric(gsub("%", "", tumor_var_freq))) %>%
  dplyr::filter(tumor_var_freq > tvf, somatic_status =="Somatic", somatic_p_value < pval)

del.tab<- 
  indel.tab %>% 
  filter(grepl("-", .$var, fixed=T)) %>%
  mutate(ref = paste0(ref, gsub("-", "", var, fixed=T)),
         var="-",
         position.end = position+nchar(ref)-1)

in.tab<- 
  indel.tab %>% 
  filter(grepl("+", .$var, fixed=T)) %>%
  mutate(var = gsub("+", "", var, fixed=T),
         ref="-",
         position=position + 1,
         position.end = position)

var.tab<- rbind(in.tab, del.tab, snp.tab) %>%
  arrange(chrom, position)

out <- var.tab[,c(1:2,24,3:4)]
fileName <- paste0(prefix, ".ann.pass")
write.table(out, file=fileName, col.names=FALSE, row.names=FALSE, quote=FALSE)

