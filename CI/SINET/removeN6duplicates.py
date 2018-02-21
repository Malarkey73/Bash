from getopt import getopt
import sys
import os

import seqUtils

HELP_STRING = """This script removes duplicates using NuGEN's N6 exra barcode technology.  
It removes reads as duplicates if they fulfill the following criteria: a) start in the same location b) have the same N6 barcode.  If reads are removed, the read with the highest mapping quality is retained.  By default the last 6 bases of the index read sequence is used.  To use a different set of bases, use the -s and -x parameters.

    -i   input SORTED SAM file for Read 1
    -b   barcode fastq file, can be the full index read or the read fastq with the barcode in the title in 1:N:0:TTTTTTTTTTTTTT format
    -t   whether or not the barcode in the title of the fastq record.  The default is False, meaning the barcode is in the read portion (the second line) of each fastq record.  To turn this on, specify "-t True"
    -o   output sorted sam file with dups removed
    -s   where to start considering the N6 (default = last 6 bases are considered)
    -x   where to stop considering the N6 (default = -1 = end of the read)
    -l   log file to write statistics to (optional)
"""

# DEFAULTS 
DEFAULT_START = -1
DEFAULT_STOP = -1
DEFAULT_LENGTH = 6
DEFAULT_IS_IN_TITLE = False
# for now, pretend all reads have the same length.  If reads have different lengths, there
# are two possible reasons -- either they were trimmed for adapter or they were trimmed for
# quality.  If they were trimmed for adapter, they are different sources, shouldn't be considered
# PCR duplicates and we should consider different lengths = different reads.  But if they were
# trimmed for quality, they could be PCR duplicates.  In this case, the latter is more likely
# so we ignore length when we determine PCR duplicates and just go on start site and N6 code.
ALL_SAME_LENGTH = 10

def writeOneLine(out, line, bc):
    #print "for bc %s, printing read %s" % (bc, line.split()[0])
    pieces = line.split("\t")
    #reads_kept.append(pieces[0])
    pieces[0] = "%s__%s" % (pieces[0], bc)
    out.write("\t".join(pieces))

def processLines(barcodes, linesToProcess, out):
    """Processes all the reads that start in one place."""
    # process the lines we've accumulated, 
    # one loop for each individual barcode
    alignmentsProcessed = 0
    newMethodDups = 0
    oldMethodDups = 0

    for (bc, allLengths) in linesToProcess.iteritems():
        for (seqlength, lines) in allLengths.iteritems():
            alignmentsProcessed += len(lines)
            if len(lines) == 1:
                #print "bc %s is a singleton" % (bc)
                writeOneLine(out, lines[0], bc)
            else:
                maxQuality = 0
                for l in lines:
                    maxQuality = max(maxQuality, int(l.split()[4]))
                #print "Max quality is %s" % (maxQuality)
                for l in lines:
                    quality = int(l.split()[4])
                    if quality == maxQuality:
                        writeOneLine(out, l, bc)
                        break
                newMethodDups += len(lines) - 1

    if alignmentsProcessed:
        oldMethodDups += alignmentsProcessed - 1
    else:
        oldMethodDups = 0

    return alignmentsProcessed, newMethodDups, oldMethodDups

# retrieve the user parameters
read1sam = None
barcodeFile = None
isInTitle = DEFAULT_IS_IN_TITLE
#inputFastq = None
outputFilename = None
startIndex = DEFAULT_START
stopIndex = DEFAULT_STOP
logFile = None

try:
    optlist, args = getopt(sys.argv[1:], "hi:b:o:s:x:l:t:")
except:
    print "Error retrieving options"
    print ""
    print HELP_STRING
    sys.exit(1)

for (opt, opt_arg) in optlist:
    if opt == "-h":
        print ""
        print HELP_STRING
        sys.exit(1)
    elif opt == "-i":
        read1sam = opt_arg
    elif opt == "-b":
        barcodeFile = opt_arg
    elif opt == "-t":
        isInTitle = opt_arg.startswith("T") or opt_arg.startswith("t")
    #elif opt == "-c":
    #    inputFastq = opt_arg
    elif opt == "-o":
        outputFilename = opt_arg
    elif opt == "-s":
        startIndex = int(opt_arg)
    elif opt == "-x":
        stopIndex = int(opt_arg)
    elif opt == "-l":
        logFile = opt_arg

# check required parameters exist
if read1sam == None or barcodeFile == None or outputFilename == None:
    print "You must specify the SORTED SAM file (-i), barcode file (-b) and output file (-o)"
    print ""
    print HELP_STRING
    sys.exit(1)
    
# do the actual work

outLog = None
if logFile:
    outLog = open(logFile, "w")

# read in the barcode indexes for each read
reads_kept = []
barcodes = {}
for (title, seq, qual) in seqUtils.FastqIterator(barcodeFile):
    #if useTitle == None:
    #    part1 = title.split()[-1]
    #    if part1.count(":") == 3 and part1.split(":")[3] != "0" and len(part1.split(":")[3]) > 0:
    #        useTitle = True
    #        seq = part1.split(":")[3]
    #        if startIndex < 0 and stopIndex < 0:
    #            bc = seq[-DEFAULT_LENGTH:]
    #        elif stopIndex < 0:
    #            bc = seq[startIndex:]
    #        else:
    #            bc = seq[startIndex:stopIndex]
    #        print """The barcode appears to be in the fastq title.  For the first record, the title is:
#%s
#The barcode is parsed from the information in the title after the first space, in this case:
#%s
#The part of the barcode that will be used for the N6 code is the characters after the final colon:
#%s
#""" % (title, part1, bc)
#        else:
#            useTitle = False
#            if startIndex < 0 and stopIndex < 0:
#                bc = seq[-DEFAULT_LENGTH:]
#            elif stopIndex < 0:
#                bc = seq[startIndex:]
#            else:
#                bc = seq[startIndex:stopIndex]
#
#            print """The fastq file appears to be the index read.  For the first record, the sequence is:
#%s
#And the barcode will be the following:
#%s
#""" % (seq, bc)

    if isInTitle:
        part1 = title.split()[-1]
        if part1.count(":") != 3:
            raise Exception("We expect the fastq title to end in something like 1:N:0:CACGTCTATAACTC, where the final piece is the barcode.  Your title ends with '%s', which doesn't fit the same pattern.")
        seq = part1.split(":")[3]
        title = title.split()[0]
        if startIndex < 0 and stopIndex < 0:
            barcodes[title] = seq[-DEFAULT_LENGTH:]
        elif stopIndex < 0:
            barcodes[title] = seq[startIndex:]
        else:
            barcodes[title] = seq[startIndex:stopIndex]
    else:
        title = title.split()[0]
        if startIndex < 0 and stopIndex < 0:
            barcodes[title] = seq[-DEFAULT_LENGTH:]
        elif stopIndex < 0:
            barcodes[title] = seq[startIndex:]
        else:
            barcodes[title] = seq[startIndex:stopIndex]

#print "Done reading barcodes.  Some examples:"
#count = 0
#for (k,v) in barcodes.iteritems():
#    print "%s = %s" % (k, v)
#    count += 1
#    if count > 10:
#        break

# read through the sorted sam file
linesToProcessFwd = {}
linesToProcessRev = {}
curChrom = None
curPosition = None

out = open(outputFilename, "w")
alignmentsProcessed = 0
newMethodDups = 0
oldMethodDups = 0
unalignedReads = 0

for line in open(read1sam):
    #if alignmentsProcessed > 500:
    #    break

    if line.startswith("@"):
        out.write(line)
        continue

    pieces = line.split("\t")
    chrom = pieces[2]
    if chrom == "*":
        unalignedReads += 1
        continue

    #alignmentsProcessed += 1
    isForward = int(pieces[8]) > 0
    position = int(pieces[3])
    bc = barcodes[pieces[0].split('_')[0]] #editted 20140815
    seqlength = ALL_SAME_LENGTH #len(pieces[9])
    if chrom == curChrom and position == curPosition:
        if isForward:
            if not linesToProcessFwd.has_key(bc):
                linesToProcessFwd[bc] = {}
            if not linesToProcessFwd[bc].has_key(seqlength):
                linesToProcessFwd[bc][seqlength] = []
            linesToProcessFwd[bc][seqlength].append(line)
        else:
            if not linesToProcessRev.has_key(bc):
                linesToProcessRev[bc] = {}
            if not linesToProcessRev[bc].has_key(seqlength):
                linesToProcessRev[bc][seqlength] = []
            linesToProcessRev[bc][seqlength].append(line)
    else:
        #print "Processing reads for %s at %s going forward:" % (curChrom, curPosition)
        #for (k, v) in linesToProcessFwd.iteritems():
        #    print "\t%s" % (k)
        #    for x in v:
        #        print "\t\t%s" % (x.split()[0])
        #print "Processing reads for %s at %s going reverse:" % (curChrom, curPosition)
        #for (k, v) in linesToProcessRev.iteritems():
        #    print "\t%s" % (k)
        #    for x in v:
        #        print "\t\t%s" % (x.split()[0])
        a, n, o = processLines(barcodes, linesToProcessFwd, out)
        alignmentsProcessed += a
        newMethodDups += n
        oldMethodDups += o
        #print "There were %s alignments forward, with %s new method dups and %s old method dups" % (a, n, o)
        a, n, o = processLines(barcodes, linesToProcessRev, out)
        alignmentsProcessed += a
        newMethodDups += n
        oldMethodDups += o
        #print "There were %s alignments reverse, with %s new method dups and %s old method dups" % (a, n, o)

        # and set up to loop
        curChrom = chrom
        curPosition = position
        linesToProcessFwd = {}
        linesToProcessRev = {}
        if isForward:
            linesToProcessFwd[bc] = {}
            linesToProcessFwd[bc][seqlength] = [line]
        else:
            linesToProcessRev[bc] = {}
            linesToProcessRev[bc][seqlength] = [line]


##process last set
a, n, o = processLines(barcodes,linesToProcessFwd,out)
alignmentsProcessed += a
newMethodDups += n
oldMethodDups += o
a, n, o = processLines(barcodes,linesToProcessRev,out)
alignmentsProcessed += a
newMethodDups += n
oldMethodDups += o

print "Done!  There were %s unaligned reads, %s alignments processed and %s duplicates (dup rate = %.2f).  Without the barcodes, there would have been %s duplicates (dup rate = %.2f)" % (unalignedReads, alignmentsProcessed, newMethodDups, ( (newMethodDups * 100.0) / alignmentsProcessed), oldMethodDups, ( (oldMethodDups * 100.0) / alignmentsProcessed) )

if outLog:
    outLog.write("Done!  There were %s unaligned reads, %s alignments processed and %s duplicates (dup rate = %.2f).  Without the barcodes, there would have been %s duplicates (dup rate = %.2f)\n" % (unalignedReads, alignmentsProcessed, newMethodDups, ( (newMethodDups * 100.0) / alignmentsProcessed), oldMethodDups, ( (oldMethodDups * 100.0) / alignmentsProcessed) ))
