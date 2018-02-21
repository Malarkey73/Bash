from types import StringTypes

def convertQualityStr(strVersion, offset=33):
    """Converts a quality string to a list of integers.                                                              
    an offset of 33 is Phred style and an offet of 64 is Solexa style"""
    return map( lambda x: ord(x)-offset, strVersion )

def convertToQualityStr(intList, offset=33):
    """Converts a list of integers to a quality string."""
    return "".join( map(lambda x: chr(x+offset), intList) )

# creates an iterator that loops through a fastq file returns records one at a time
def FastqIterator(fh):
    """return an iterator of Records found in file handle, fh.
    """
    def readTotitle(fh, titleChar):
        """returns a tuple ([lines before the next title line], next tile line)
        """
        preLines = []
        while True:
            l = fh.readline().strip()
            if l.startswith(titleChar):
                return (preLines,l)
            elif l == '':
                return preLines,None
            else:
                preLines.append(l)

    if type(fh) in StringTypes:
        fh = file(fh)
    
    preLines,nextTitleLine =readTotitle(fh,'@')

    while nextTitleLine != None:
        seqTitle = nextTitleLine[1:].rstrip()
        preLines,nextTitleLine=readTotitle(fh,'+')
        qualTitle = nextTitleLine[1:].rstrip()
        if len(qualTitle.strip()) > 0 and seqTitle != qualTitle:
            print seqTitle
            print preLines
            print qualTitle
            raise Exception("Error in parsing: @title sequence entry must be immediately followed by corresponding +title quality entry.")
        seqLines = preLines
        qualLines = []
        for i in range(len(seqLines)): # Quality characters should be the same length as the sequence
            qualLines.append( fh.readline().strip() )
        
        preLines,nextTitleLine=readTotitle(fh,'@')

        yield (seqTitle, ''.join(seqLines), ''.join(qualLines))


def FastaIterator(fh):
    """return an iterator of Records found in file handle, fh.
        Original code from Kael and Dale.
    """
    def readTotitle(fh):
        """returns a tuple ([lines before the next title line], next tile line)
        """
        preLines = []
        while True:
            l = fh.readline().strip()
            if l.startswith('>'):
                return (preLines,l)
            elif l == '':
                return preLines,None
            else:
                preLines.append(l)


    if type(fh) in StringTypes:
        fh = file(fh)
    
    preLines,nextTitleLine =readTotitle(fh)

    while nextTitleLine != None:
        title = nextTitleLine[1:].rstrip()
        preLines,nextTitleLine=readTotitle(fh)
        yield (title,''.join(map(lambda x: x.rstrip(),preLines)))


def fastqToFasta(fastqFile, fastaFile):
    out = open(fastaFile, "w")

    numRec = 0
    for (title, sequence, quality) in FastqIterator(fastqFile):
        out.write(">%s\n%s\n" % (title, sequence))
        numRec += 1

    print "Done.  Converted %s records." % (numRec)


def trimFastq(inputFastq, outputFastq, trimFront, trimRear):
    """Creates a new fastq file with the given number of bases removed from the
    front or rear of the read."""

    if trimFront == 0 and trimRear == 0:
        raise Exception("trimFront and trimRear are both zero.  Nothing to trim.")

    out = open(outputFastq, "w")
    numRec = 0
    for (title, sequence, quality) in FastqIterator(inputFastq):
        if numRec == 0 and (trimFront + trimRear) > len(sequence):
            raise Exception("The sequence is shorter than the trim lengths!")
        if trimFront > 0 and trimRear > 0:
            subseq = sequence[trimFront:-trimRear]
            subqual = quality[trimFront:-trimRear]
        elif trimFront > 0:
            subseq = sequence[trimFront:]
            subqual = quality[trimFront:]
        else:
            subseq = sequence[:-trimRear]
            subqual = quality[:-trimRear]
        #print "%s\t%s\t%s" % (title, subseq, subqual)
        #print len(sequence), len(subseq)
        out.write("@%s\n%s\n+\n%s\n" % (title, subseq, subqual))
        numRec += 1
        #if numRec > 10:
        #    break

    print "Done.  There were %s records that were trimmed." % (numRec)


def fastaToDict(inputFasta):

    d = {}
    for (title, sequence) in FastaIterator(inputFasta):
        d[title] = sequence

    return d



def makeGenomeTxtFile(inputFasta, outputGenomeFile):
    
    out = open(outputGenomeFile, "w")

    for (title, sequence) in FastaIterator(inputFasta):
        out.write("%s\t%s\n" % (title, len(sequence)))


def calculateGCpercentForSeq(seq):
    
    totalAT =0 
    totalGC = 0
    for base in seq:
        if base == "A" or base == "T" or base == "a" or base == "t":
            totalAT += 1
        elif base == "C" or base == "G" or base == "c" or base == "g":
            totalGC += 1
        elif base == "N" or base == "n":
            pass
        else:
            raise Exception("Illegal base '%s' found in '%s'" % (base, seq))

    # hmmm.... what's the right thing to do for bins that are all N?
    if totalAT+totalGC == 0:
        return 0
    return ( (totalGC * 100.0) / (totalAT+totalGC))

def calculateGCpercentForFastq(inputFastq):

    totalAT = 0
    totalGC = 0
    otherBreakdown = {}
    for (title, sequence, quality) in FastqIterator(inputFastq):
        for base in sequence:
            if base == "A" or base == "T":
                totalAT += 1
            elif base == "C" or base == "G":
                totalGC += 1
            else:
                #print "ACK!  Non A,C,G,T base detected: %s" % (base)
                if not otherBreakdown.has_key(base):
                    otherBreakdown[base] = 0
                otherBreakdown[base] += 1

    print "There were %s AT bases, %s GC bases and %s other bases" % (totalAT, totalGC, sum(otherBreakdown.values()))
    print "Other bases were:"
    for (k,v) in otherBreakdown.iteritems():
        print "\t%s\t%s" % (k, v)
    print "The percent GC = %.3f   (with the other bases not considered!)" % ( totalGC / float(totalAT + totalGC))


def readTestFastq(inputFastq, printEvery=100000, offset=33, breakOnQualError=True):
    """Simply reads through a fastq file, printing the current line number every once
    in a while, to try to detect corrupted lines in the file."""

    count = 0
    countBad = 0
    for (title, sequence, quality) in FastqIterator(inputFastq):
        count += 1
        if count % printEvery == 0:
            print count, title

        qualInts = convertQualityStr(quality, offset)
        for i in range(len(qualInts)):
            q = qualInts[i]
            if q < -5 or q > 45:
                print "The quality appears to be wrong in record '%s'. (count=%s)" % (title, count)
                print "In the quality string '%s', the value '%s' converts to quality score %s" % (quality, quality[i], q)
                countBad += 1
                if breakOnQualError:
                    raise Exception("The quality appears to be wrong in record '%s'.\nIn the quality string '%s', the value '%s' converts to quality score %s"% (title, quality, quality[i],q))

        if sequence.find(".") >= 0:
            print "The record '%s' has sequence '%s' which has a period in it.  Convert these to Ns to continue. (count=%s)" % (title, sequence, count)
            countBad += 1
            if breakOnQualError:
                raise Exception("The record '%s' has sequence '%s' which has a period in it.  Convert these to Ns to continue." % (title, sequence))
        
        if len(sequence) != len(quality):
            print "The record '%s' has a sequence of length %s but a quality of length %s (count=%s)" % (title, len(sequence), len(quality), count)
            countBad += 1
            if breakOnQualError:
                raise Exception("The record '%s' has a sequence of length %s but a quality of length %s" % (title, len(sequence), len(quality)))

    print "Done!  There were %s records and %s bad records" % (count, countBad)

def reverseComplement(seq):
    
    seq = seq[::-1]
    retval = []
    for base in seq:
        if base == "A":
            retval.append("T")
        elif base == "C":
            retval.append("G")
        elif base == "G":
            retval.append("C")
        elif base == "T":
            retval.append("A")
        elif base == "a":
            retval.append("t")
        elif base == "c":
            retval.append("g")
        elif base == "g":
            retval.append("c")
        elif base == "t":
            retval.append("a")
        elif base == "N":
            retval.append("N")
        else:
            raise Exception("Illegal base %s in reverseComplement of %s" % (base, seq))

    return "".join(retval)


def fixFasta(inFasta, outFasta, lineLength=60):
    """Given a fasta file with uneven sequence lines, fixes the sequence lines so they
    are all lengthLength (default 60) long."""

    out = open(outFasta, "w")
    for (title, sequence) in FastaIterator(inFasta):
        out.write(">%s\n" % (title))
        for i in range(0, len(sequence), lineLength):
            out.write("%s\n" % sequence[i:i+lineLength])
        

