require(ggplot2)
theme(theme_bw())
# the names od each type of files
ccurveFiles=dir()[grep("ccurve.txt", dir())]
lcextrapFiles=dir()[grep("lcextrap.txt", dir())]

for(ccF in ccurveFiles){
  ccurve=read.table(cc, header=T)
  qplot(TOTAL_READS, DISTINCT_READS, data=ccurve)
  ggsave(gsub(".txt", ".pdf", ccF))
}

for(lcexF in lcextrapFiles){
  lcextrap=read.table(lcexF, header=T)
  qplot(TOTAL_READS, EXPECTED_DISTINCT, data=lcextrap)+
    geom_ribbon(aes(ymin=LOGNORMAL_LOWER_95.CI, ymax=LOGNORMAL_UPPER_95.CI), alpha=.2)
  ggsave(gsub(".txt", ".pdf"))
}
