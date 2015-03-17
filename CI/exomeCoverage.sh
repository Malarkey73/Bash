# Set up to run in GNU parallel.
find *bam | parallel 'bedtools coverage -hist -abam {} -b target_regions.bed | grep ^all > {}.hist.all.txt'
 
# What this actually runs:
bedtools coverage -hist -abam samp.01.bam -b target_regions.bed | grep ^all > samp.01.bam.hist.all.txt
bedtools coverage -hist -abam samp.02.bam -b target_regions.bed | grep ^all > samp.02.bam.hist.all.txt
bedtools coverage -hist -abam samp.03.bam -b target_regions.bed | grep ^all > samp.03.bam.hist.all.txt
bedtools coverage -hist -abam samp.05.bam -b target_regions.bed | grep ^all > samp.05.bam.hist.all.txt
bedtools coverage -hist -abam samp.06.bam -b target_regions.bed | grep ^all > samp.06.bam.hist.all.txt
bedtools coverage -hist -abam samp.10.bam -b target_regions.bed | grep ^all > samp.10.bam.hist.all.txt