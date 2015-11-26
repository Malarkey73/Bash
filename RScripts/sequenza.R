#!/usr/bin/Rscript

# usage 
# here the assumption is existing matched NB7.snp and NB7.copynumber files
# nohup ./sequenza.R NB7

args<-commandArgs(trailingOnly=TRUE)
prefix<-args[1]
#read and convert Varscan output
snp <- read.table(paste0(prefix,".snp"), header = TRUE, sep = "\t")
cnv <- read.table(paste0(prefix,".copynumber"), header = TRUE, sep = "\t")

library(sequenza)
seqz.data <- VarScan2seqz(varscan.somatic = snp, varscan.copynumber = cnv)

# gc normalisation
gc.stats <- gc.norm(x = seqz.data$depth.ratio, gc = seqz.data$GC.percent)
gc.vect  <- setNames(gc.stats$raw.mean, gc.stats$gc.values)
gc.smoothed<- with(gc.stats, loess(raw.mean~gc.values))
gc.vect <-setNames(gc.smoothed$fitted, gc.stats$gc.values)
seqz.data$adjusted.ratio <- seqz.data$depth.ratio/gc.vect[as.character(seqz.data$GC.percent)]

# plot normalisation
pdf(file=paste0(prefix, "_gcnorm.pdf"))
par(mfrow = c(1,2), cex = 1, las = 1, bty = 'l')
with(gc.stats, plot(gc.values, raw.mean))
with(gc.smoothed, lines(x[,1], fitted, col="red", lwd=3))
hist2(seqz.data$depth.ratio, seqz.data$adjusted.ratio, 
      breaks = prettyLog, key = vkey, panel.first = abline(0, 1, lty = 2),
      xlab = 'Uncorrected depth ratio', ylab = 'GC-adjusted depth ratio')
dev.off()

# seqz fannying about
out <- paste0(prefix, ".seq.gz")
write.table(seqz.data, gzfile(out), col.names = TRUE, row.names = FALSE, sep = "\t")
extract <- sequenza.extract(out)

# get the dat awe actually need out
er<-do.call(rbind, extract$ratio)
er$chr<- dirname(gsub(".", "/", row.names(er), fixed=T))
seg<-do.call(rbind, extract$segments)
write.table(er, file=paste0(prefix, "_ratios.txt"), quote=F, sep="\t")
write.table(seg, file=paste0(prefix, "_segments.txt"), quote=F, sep="\t")

# hint for how to handle that output
#  er <- read.delim("NB1_ratios.txt, header=T)
#  seg<- read.delim("NB1_segments.txt, header=T)
#  qplot(x=start, xend=end, y=mean, yend=mean, data=er, geom="segment")+
#  theme_bw()+ 
#  geom_ribbon(aes(ymin=q0, ymax=q1), alpha=0.1)+
#  geom_segment(aes(x=start.pos, xend=end.pos, y=depth.ratio, yend=depth.ratio), data=seg, col="red")




