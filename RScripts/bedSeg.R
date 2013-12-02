  args <- commandArgs(trailingOnly = TRUE)
  bedGraph=as.character(args[1])
  segWidth=as.numeric(args[2])
  require(dplyr)
  require(data.table)
  # out will be name of file output
  out=paste0(bedGraph,'.seg')
  # to speed up long read i use bash to get the number of rows so mem can be preallocated
  numRow=as.integer(system(paste("wc -l", bedGraph, "| sed 's/[^0-9.]*\\([0-9.]*\\).*/\\1/'"), intern=T))
  # data.table fread is blazingly fast
  BG=fread(bedGraph, nrows=numRow, header=F, colClasses=c('character', rep('integer',3)))
  setnames(BG, colnames(BG), c('chr', 'start', 'end', 'score'))
  # zikes plyr smartness adds bin colummns - one pass through data
  BG=mutate(BG, startSeg=floor(((start+end)/2)/segWidth)*segWidth, endSeg=startSeg+500)
  # and aggregate for per chr per bin score
  BG=summarise(group_by(BG, chr, startSeg, endSeg), score=mean(score))
  write.table(BG, file=out, sep='\t', quote=F, row.names=F, col.names=F)
