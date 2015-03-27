#!/usr/bin
ANNOTATESCRIPT="/home/rmgzshd/annovar/annotate_variation.pl"
ANNOTATEDB="/home/rmgzshd/annovar/humandb"

export ANNOTATESCRIPT; export ANNOTATEDB

# this takes all the .pass files and converts them to .var annovar input files
my_awk='{print $1,$2,$2,$3,$4}'
parallel "cat {} | awk '$my_awk' > {.}.var" ::: *.pass

# this runs annovar on them all at once
ls *.var | parallel $ANNOTATE -out {.} -build hg19 {} $ANNOTATEDB