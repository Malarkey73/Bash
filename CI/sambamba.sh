#!/bin/bash

# this is an 18GB file

time(sambamba slice TCGA-CN-5367-01A-01D-1431_120420_SN1120_0134_AC0J8YACXX_s_2_rg.sorted.bam 1:100000-110000 |
sambamba view -F "unmapped or mate_is_unmapped" /dev/stdin > test1.bam)
#real	0m0.173s
#user	0m0.151s
#sys	0m0.205s

time(sambamba view -F "unmapped or mate_is_unmapped" TCGA-CN-5367-01A-01D-1431_120420_SN1120_0134_AC0J8YACXX_s_2_rg.sorted.bam 1:100000-110000 > test1.bam)
#real	0m0.196s
#user	0m0.227s
#sys	0m0.072s

time(sambamba slice TCGA-CN-5367-01A-01D-1431_120420_SN1120_0134_AC0J8YACXX_s_2_rg.sorted.bam 15:100000-110000 |
sambamba view -F "unmapped or mate_is_unmapped" /dev/stdin > test3.bam)
#real	0m0.145s
#user	0m0.133s
#sys	0m0.096s

time(sambamba view -F "unmapped or mate_is_unmapped" TCGA-CN-5367-01A-01D-1431_120420_SN1120_0134_AC0J8YACXX_s_2_rg.sorted.bam 15:100000-110000 > test4.bam)
#real	0m0.224s
#user	0m0.194s
#sys	0m0.115s


time(sambamba view -F "unmapped or mate_is_unmapped" TCGA-CN-5367-01A-01D-1431_120420_SN1120_0134_AC0J8YACXX_s_2_rg.sorted.bam > test5.bam)
#real	1m13.549s
#user	12m3.303s
#sys	0m29.538s



time(
	sambamba view  -F  "unmapped or mate_is_unmapped" TCGA-CN-5367-01A-01D-1431_120420_SN1120_0134_AC0J8YACXX_s_2_rg.sorted.bam |
	grep -v ^@ | 
	tee >(awk 'NR%2==1 {print "@"$1"\n"$10"\n+\n"$11}' > test.fq1) > >(awk 'NR%2==0 {print "@"$1"\n"$10"\n+\n"$11}' > test.fq2)
	)
#real	2m54.136s
#user	17m40.941s
#sys	1m1.604s


time(
	$RETROSEQ -discover -bam $TARGETBAM -output $prefix.candidates.tab -refTEs $REFS -align
	)

time(
	$RETROSEQ -discover -bam $TARGETBAM -output $prefix.candidates.tab -eref $PROBES -align
	)




