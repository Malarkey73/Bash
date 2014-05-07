#!/bin/bash
SAMTOOLS="/home/rmgzshd/samtools/samtools"
$SAMTOOLS merge mus_ChIP.bam ERR022291.bam ERR022305.bam
$SAMTOOLS merge mus_Inp.bam ERR022272.bam ERR022273.bam ERR022296.bam
$SAMTOOLS merge TC1_ChIP.bam ERR022297.bam ERR022298.bam ERR022299.bam
$SAMTOOLS merge TC1_Inp.bam ERR022278.bam ERR022279.bam ERR022280.bam ERR022281.bam
