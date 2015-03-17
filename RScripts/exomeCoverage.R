# Assumes you've already run coverageBed -hist, and grep'd '^all'. E.g. something like:
# find *.bam | parallel 'bedtools -abam {} -b capture.bed -hist | grep ^all > {}.all.txt'
 
# Get a list of the bedtools output files you'd like to read in
print(files <- list.files(pattern="all.txt$"))
 
# Optional, create short sample names from the filenames. 
# For example, in this experiment, my sample filenames might look like this:
# prefixToTrash-01.pe.on.pos.dedup.realigned.recalibrated.bam
# prefixToTrash-02.pe.on.pos.dedup.realigned.recalibrated.bam
# prefixToTrash-03.pe.on.pos.dedup.realigned.recalibrated.bam
# This regular expression leaves me with "samp01", "samp02", and "samp03" in the legend.
print(labs <- paste("samp", gsub("prefixToTrash-0|\\.pe\\.on\\.pos\\.dedup\\.realigned\\.recalibrated\\.bam\\.cov\\.hist\\.txt\\.all\\.txt", "", files, perl=TRUE), sep=""))
 
# Create lists to hold coverage and cumulative coverage for each alignment,
# and read the data into these lists.
cov <- list()
cov_cumul <- list()
for (i in 1:length(files)) {
    cov[[i]] <- read.table(files[i])
    cov_cumul[[i]] <- 1-cumsum(cov[[i]][,5])
}
 
# Pick some colors
# Ugly:
# cols <- 1:length(cov)
# Prettier:
# ?colorRampPalette
# display.brewer.all()
library(RColorBrewer)
cols <- brewer.pal(length(cov), "Dark2")
 
# Save the graph to a file
png("exome-coverage-plots.png", h=1000, w=1000, pointsize=20)
 
# Create plot area, but do not plot anything. Add gridlines and axis labels.
plot(cov[[1]][2:401, 2], cov_cumul[[1]][1:400], type='n', xlab="Depth", ylab="Fraction of capture target bases \u2265 depth", ylim=c(0,1.0), main="Target Region Coverage")
abline(v = 20, col = "gray60")
abline(v = 50, col = "gray60")
abline(v = 80, col = "gray60")
abline(v = 100, col = "gray60")
abline(h = 0.50, col = "gray60")
abline(h = 0.90, col = "gray60")
axis(1, at=c(20,50,80), labels=c(20,50,80))
axis(2, at=c(0.90), labels=c(0.90))
axis(2, at=c(0.50), labels=c(0.50))
 
# Actually plot the data for each of the alignments (stored in the lists).
for (i in 1:length(cov)) points(cov[[i]][2:401, 2], cov_cumul[[i]][1:400], type='l', lwd=3, col=cols[i])
 
# Add a legend using the nice sample labeles rather than the full filenames.
legend("topright", legend=labs, col=cols, lty=1, lwd=4)
 
dev.off()