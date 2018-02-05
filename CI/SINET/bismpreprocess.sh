#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

TRIMGALORE="/home/rmgzshd/TrimGalore-0.4.3/trim_galore"
TRIMDIVERSITY="/mnt/store1/TESTBISMARK/trimRRBSdiversityAdaptCustomers.py"
PYTHON27="/usr/local/bin/python2.7"
NUDUP="/mnt/store1/TESTBISMARK/nudup.py"



while read COL1
	do

	iget $COL1
	$TRIMGALORE --dont_gzip -a AGATCGGAAGAGC $COL1





	done < GETFQFROMIRODS