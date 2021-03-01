#!/usr/bin/env python

import sys
import gzip
import os.path
import Bio.pairwise2 as pw2

PRIMER1 = "CGAGCGGGTC"
PRIMER2 = "CGGAGGTCAT"

SGRNA1 = "CACCACCGCCAAGATCACCA"

#                *          * 
SGRNA2   = "GGTCANTCTGCACAGCNCCA"
SGRNA2v1 = "GGTCAATCTGCACAGCGCCA"  # F F 0
SGRNA2v2 = "GGTCAATCTGCACAGCACCA"  # F T 1
SGRNA2v4 = "GGTCAGTCTGCACAGCGCCA"  # T F 2
SGRNA2v3 = "GGTCAGTCTGCACAGCACCA"  # T T 3

# Conserved sequences
CSEQ1 = "ACCGCGTCTGC"
CSEQ2 = "TGAAGCTTCTC"

SUB = 1
INS = 2
DEL = 4
LONG = 8
SHORT = 16
COUNTERS = [ 0,   1,   2,   3,    4,   5,    6,    7,     14,    15,     22,     23]
COLNAMES = ["M", "S", "I", "SI", "D", "SD", "ID", "SID", "ID*", "SID*", "ID**", "SID**"]

""" Possible values:
0 M
1 S
2 I
3 SI
4 D
5 SD
6 ID
7 SID

14 ID*
15 SID*

22 ID**
23 SID**

30 IDN
31 SIDN
"""

def averageQual(qstring):
    n = 0
    for c in qstring:
        n += ord(c) - 33
    return 1.0*n/len(qstring)

class Read(object):
    seq = ""
    qual = ""
    avgq = 0

    def __init__(self, seq, qual):
        self.seq = seq
        self.qual = qual
        self.avgq = averageQual(qual)

class Averager(object):
    s = 0
    n = 0

    def __init__(self):
        self.s = 0
        self.n = 0
    
    def add(self, x):
        self.s += x
        self.n += 1

    def avg(self):
        if self.n == 0:
            return 0.0
        else:
            return 1.0*self.s/self.n

class Counter(object):
    target = ""
    tglen = 0
    nhits = 0
    counters = []
    variants = []

    cs_good = 0
    cs_ins = 0
    cs_ins_avg = None
    cs_del = 0
    cs_del_avg = None
    
    debug = False
    
    def __init__(self, target):
        self.target = target
        self.tglen = len(target)
        self.counters = [0]*24
        self.variants = [ {'A': 0, 'C': 0, 'G': 0, 'T':0, 'N':0} for i in range(self.tglen) ]
        self.inslen = {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0}
        self.dellen = {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0}
        self.insavg = Averager()
        self.delavg = Averager()

        self.cs_good = 0
        self.cs_ins = 0
        self.cs_ins_avg = Averager()
        self.cs_del = 0
        self.cs_del_avg = Averager()
        
    def checkTarget(self, read):
        hits = pw2.align.localms(read, self.target, 1, -.5, -.5, -.2)
        if hits:
            (seq, tgseq, score, start, end) = hits[0]
            if score > 15:
                self.update(seq, tgseq, start, end)
                return True

        # No-hit read
        return False

    def checkConsSeq(self, read, cs, pos):
        hits = pw2.align.localms(read, cs, 1, -1, -2, -2)
        if not hits:
            return
        (hread, hcs, score, p, end) = hits[0]
        if score < 9:
            return

        self.cs_good += 1
        if p > pos:
            self.cs_ins += 1
            #print "I" + str(p-pos)
            self.cs_ins_avg.add(p-pos)
        elif p < pos:
            self.cs_del += 1
            #print "D" + str(pos-p)
            self.cs_del_avg.add(pos-p)
    
    def update(self, seq, tgseq, start, end):
        if seq[start:start+2] != self.target[0:2]:
            return
        self.nhits += 1
        hitlen = end - start
        # Fix end if last base is a mismatch
        if hitlen == self.tglen - 1 and seq[end] != "-" and tgseq[end] != "-":
            end += 1
            hitlen += 1

        idx = 0
        fseq = seq[start:end]
        ftar = tgseq[start:end]
        ndel = fseq.count("-")  # number of deletions
        nins = ftar.count("-")  # number of insertions
        delta = nins - ndel
        cond = (nins!=0)*2+(ndel!=0) # 0=match, 1=del only, 2=ins only, 3=both

        if cond == 1:
            idx += DEL
            self.delavg.add(ndel)
            ndel = min(ndel, 6)
            self.dellen[ndel] += 1
        elif cond == 2:
            idx += INS
            self.insavg.add(nins)
            nins = min(nins, 6)
            self.inslen[nins] += 1
        elif cond == 3:
            idx += INS+DEL
            if delta > 0:
                idx += LONG
                self.insavg.add(delta)
                nins = min(delta, 6)
                self.inslen[nins] += 1
            elif delta < 0:
                idx += SHORT
                self.delavg.add(-delta)
                ndel = min(-delta, 6)
                self.dellen[ndel] += 1

        if self.hasSubst(fseq, ftar):
            idx += SUB

        self.counters[idx] += 1
#        print "{}\n{}\n{}\n".format(seq, tgseq, idx)
#        raw_input()
        
    def hasSubst(self, a, b):
        # a is read
        # b is target
        if self.debug:
            sys.stdout.write(a + "\n" + b + "\n")
        subst = False
        p = 0
        for i in range(len(a)):
            if a[i] == '-' or b[i] == '-':
                pass
            elif a[i] != b[i]:
                subst = True
                x = a[i]
                self.variants[p][x] += 1
                if self.debug:
                    sys.stdout.write("{}={} ".format(p, x))
            if b[i] != '-':
                p += 1
        if self.debug:
            sys.stdout.write("\n\n")
        #    raw_input()
        return subst

    def writeMatrix(self, out):
        out.write("Ref\tA\tC\tG\tT\tTot\tA%\tC%\tG%\tT%\n")
        r = 2
        for i in range(self.tglen):
            row = self.variants[i]
            out.write("{}\t{}\t{}\t{}\t{}".format(self.target[i], row['A'], row['C'], row['G'], row['T']))
            out.write("\t=sum(B{}:E{})".format(r, r))
            if row['A'] + row['C'] + row['G'] + row['T']:
                out.write("\t=B{}/F{}\t=C{}/F{}\t=D{}/F{}\t=E{}/F{}\n".format(r, r, r, r, r, r, r, r))
            else:
                out.write("\t\t\t\t\n")
            r += 1

class Counter2(Counter):
    
    def __init__(self, target):
        self.target = target
        self.tglen = len(target)
        self.nhits = [0]*4
        self.counters = [ [0]*24 for _ in range(4) ]
        self.inslen = [ {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0} for _ in range(self.tglen) ]
        self.dellen = [ {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0} for _ in range(self.tglen) ]
        self.insavg = [ Averager() for _ in range(self.tglen) ]
        self.delavg = [ Averager() for _ in range(self.tglen) ]
        #self.variants = [ {'A': 0, 'C': 0, 'G': 0, 'T':0} for i in range(self.tglen) ]

        self.cs_good = 0
        self.cs_ins = 0
        self.cs_ins_avg = Averager()
        self.cs_del = 0
        self.cs_del_avg = Averager()

    def checkTarget(self, read):
        hits = pw2.align.localms(read, self.target, 1, -.5, -.5, -.2)
        if hits:
            (seq, tgseq, score, start, end) = hits[0]
            if score > 15:
                self.update(seq, tgseq, start, end)
                #if idx:
                #print (seq[start:end], idx)
                #raw_input()
                return True
        return False

    def hasSubst(self, a, b, pv1, pv2):
        for i in range(len(a)):
            if a[i] == '-' or b[i] == '-':
                continue
            if i == pv1 or i == pv2:
                continue
            if a[i] != b[i]:
                return True
        return False
                
    def update(self, seq, tgseq, start, end):
        if seq[start:start+2] != self.target[0:2]:
            return
        
        hitlen = end - start
        # Fix end if last base is a mismatch
        if hitlen == self.tglen - 1 and seq[end] != "-" and tgseq[end] != "-":
            end += 1
            hitlen += 1
            
        # now find out which variant we're seeing
        matchseq = seq[start:end]
        var = 0
        if "AGT" in matchseq:
            var += 2
            pv1 = matchseq.find("AGT") + 1
        else:
            pv1 = matchseq.find("AAT") + 1

        if "CACC" in matchseq:
            var += 1
            pv2 = matchseq.find("CACC") + 1
        else:
            pv2 = matchseq.find("CGCC") + 1

        self.nhits[var] += 1

        idx = 0
        fseq = seq[start:end]
        ftar = tgseq[start:end]
        ndel = fseq.count("-")  # number of deletions
        nins = ftar.count("-")  # number of insertions
        delta = nins - ndel
        cond = (nins!=0)*2+(ndel!=0) # 0=match, 1=del only, 2=ins only, 3=both

        if cond == 1:
            idx += DEL
            self.delavg[var].add(ndel)
            ndel = min(ndel, 6)
            self.dellen[var][ndel] += 1
        elif cond == 2:
            idx += INS
            self.insavg[var].add(nins)
            nins = min(nins, 6)
            self.inslen[var][nins] += 1
        elif cond == 3:
            idx += INS+DEL
            if delta > 0:
                idx += LONG
                self.insavg[var].add(delta)
                nins = min(delta, 6)
                self.inslen[var][nins] += 1
            elif delta < 0:
                idx += SHORT
                self.delavg[var].add(-delta)
                ndel = min(-delta, 6)
                self.dellen[var][ndel] += 1

        if self.hasSubst(fseq, ftar, pv1, pv2):
            idx += SUB

        #print (matchseq, var, pv1, pv2, idx)
        #raw_input()
            
        self.counters[var][idx] += 1
            
class Scanner(object):
    filename = ""
    ctr1 = None
    ctr2 = []
    nreads = 0
    nbad = 0
    header = True
    averagers = {}
    _open = open

    def __init__(self, filename):
        self.filename = filename
        self.ctr1 = Counter(SGRNA1)
        self.ctr2 = Counter2(SGRNA2)
        #self.ctr2 = [Counter(SGRNA2), Counter(SGRNA2v2), Counter(SGRNA2v4), Counter(SGRNA2v3)] # order is correct!
        #self.ctr2 = Counter(SGRNA2)
        self.averagers = {'h1': Averager(), 'nh1': Averager(), 'h2': Averager(), 'nh2': Averager()}
        if os.path.splitext(filename)[1] == '.gz':
            self._open = gzip.open

    def scanFastq(self):
        ln = 0
        with self._open(self.filename, "r") as f:
            while True:
                line = f.readline()
                if not line:
                    break
                seq = f.readline().rstrip('\r\n')
                f.readline()
                qual = f.readline().rstrip('\r\n')
                r = Read(seq, qual)
                self.scanReadCons(r)

    def scanReadCons(self, read, checkCons=True):
        self.nreads += 1
        line = read.seq
        p1 = line.find(PRIMER1)
        if p1 > 0 and p1 + 65 < len(line):
            if self.ctr1.checkTarget(line[p1+20:p1+65]):
                self.averagers['h1'].add(read.avgq)
            else:
                self.averagers['nh1'].add(read.avgq)
                if checkCons:
                    self.ctr1.checkConsSeq(line[p1:], CSEQ1, 56)
                
        else:
            p2 = line.find(PRIMER2)
            if p2 > 0 and p2 + 45 < len(line):
                if self.ctr2.checkTarget(line[p2+5:p2+50]):
                    self.averagers['h2'].add(read.avgq)
                else:
                    self.averagers['nh2'].add(read.avgq)
                    if checkCons:
                        self.ctr2.checkConsSeq(line[p2:], CSEQ2, 61)

    def scanRead(self, line):
        self.nreads += 1
        p1 = line.find(PRIMER1)
        if p1 > 0 and p1 + 65 < len(line):
            self.ctr1.checkTarget(line[p1+20:p1+65])
        else:
            p2 = line.find(PRIMER2)
            if p2 > 0 and p2 + 45 < len(line):
                self.ctr2.checkTarget(line[p2+5:p2+50])

    def run0(self, header):
        self.scanFastq()
        self.writeResults()
        self.writeMatrices()

    def run(self, header):
        self.scanFastq()
        self.writeResults6(header)

    def run3(self, m):
        self.scanFastq()
        m.results[2].append(self.ctr1)
        m.results[3].append(self.ctr2)

    def writeResults(self):
        """Default results: target 1, target2."""
        if self.header:
            sys.stdout.write("Filename\tReads\tHit1\tM1\tS1\tI1\tSI1\tD1\tSD1\tID1\tSID1\tHit2\tM2\tS2\tI2\tSI2\tD2\tSD2\tID2\tSID2\n")
            self.header = False
        sys.stdout.write("{}\t{}\t{}".format(self.filename, self.nreads, self.ctr1.nhits))
        for c in self.ctr1.counters:
            sys.stdout.write("\t" + str(c))
        sys.stdout.write("\t" + str(self.ctr2.nhits))
        for c in self.ctr2.counters:
            sys.stdout.write("\t" + str(c))
        sys.stdout.write("\n")

    def writeResults2(self, header):
        """Default results for four variants of target2."""
        if header:
            sys.stdout.write("Filename\tReads\tHit1\tM1\tS1\tI1\tSI1\tD1\tSD1\tID1\tSID1\tHit2\tM2\tS2\tI2\tSI2\tD2\tSD2\tID2\tSID2\tHit3\tM3\tS3\tI3\tSI3\tD3\tSD3\tID3\tSID3\tHit4\tM4\tS4\tI4\tSI4\tD4\tSD4\tID4\tSID4\n")
        sys.stdout.write("{}\t{}".format(self.filename, self.nreads))
        for i in [0, 1, 3, 2]:
            sys.stdout.write("\t" + str(self.ctr2.nhits[i]))
            for c in self.ctr2.counters[i]:
                sys.stdout.write("\t" + str(c))
        sys.stdout.write("\n")

    def writeResults3(self, header):
        """Results for insertion lengths."""
        if header:
            sys.stdout.write("Filename\tReads\tHits1\tAI1\tI1*1\tI1*2\tI1*3\tI1*4\tI1*5\tI1+5\tHits2\tAI2\tI2*1\tI2*2\tI2*3\tI2*4\tI2*5\tI2+5\tHits3\tAI3\tI3*1\tI3*2\tI3*3\tI3*4\tI3*5\tI3+5\tHits4\tAI4\tI4*1\tI4*2\tI4*3\tI4*4\tI4*5\tI4+5\n")
        sys.stdout.write("{}\t{}".format(self.filename, self.nreads))
        for i in [0, 1, 3, 2]:
            sys.stdout.write("\t{}\t{}".format(self.ctr2.nhits[i], self.ctr2.insavg[i].avg()))
            for z in [1, 2, 3, 4, 5, 6]:
                sys.stdout.write("\t" + str(self.ctr2.inslen[i][z]))
        sys.stdout.write("\n")
            
    def writeResults4(self, header):
        """Results for deletion lengths."""
        if header:
            sys.stdout.write("Filename\tReads\tHits1\tAD1\tD1*1\tD1*2\tD1*3\tD1*4\tD1*5\tD1+5\tHits2\tAD2\tD2*1\tD2*2\tD2*3\tD2*4\tD2*5\tD2+5\tHits3\tAD3\tD3*1\tD3*2\tD3*3\tD3*4\tD3*5\tD3+5\tHits4\tAD4\tD4*1\tD4*2\tD4*3\tD4*4\tD4*5\tD4+5\n")
        sys.stdout.write("{}\t{}".format(self.filename, self.nreads))
        for i in [0, 1, 3, 2]:
            sys.stdout.write("\t{}\t{}".format(self.ctr2.nhits[i], self.ctr2.delavg[i].avg()))
            for z in [1, 2, 3, 4, 5, 6]:
                sys.stdout.write("\t" + str(self.ctr2.dellen[i][z]))
        sys.stdout.write("\n")

    def writeResults5(self, header):
        """Results for insertion and deletion lengths, sgRNA1."""
        if header:
            sys.stdout.write("Filename\tReads\tHits1\tAI1\tI1*1\tI1*2\tI1*3\tI1*4\tI1*5\tI1+5\tAD1\tD1*1\tD1*2\tD1*3\tD1*4\tD1*5\tD1+5\n")
        sys.stdout.write("{}\t{}".format(self.filename, self.nreads))
        sys.stdout.write("\t{}\t{}".format(self.ctr1.nhits, self.ctr1.insavg.avg()))
        for z in [1, 2, 3, 4, 5, 6]:
            sys.stdout.write("\t" + str(self.ctr1.inslen[z]))
        sys.stdout.write("\t{}\t{}".format(self.ctr1.nhits, self.ctr1.delavg.avg()))
        for z in [1, 2, 3, 4, 5, 6]:
            sys.stdout.write("\t" + str(self.ctr1.dellen[z]))
        sys.stdout.write("\n")

    def writeResults6(self, header):
        if header:
            sys.stdout.write("Filename\tReads\tNoHits1\tCons1\tExpected1\tIns1\tAvgIns1\tDel1\tAvgDel1\tNoHits2\tCons\tExpected2\tIns2\tAvgIns2\tDel2\tAvgDel2\n")
        sys.stdout.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
            self.filename, self.nreads,
            self.ctr1.nhits, self.ctr1.cs_good, self.ctr1.cs_good - self.ctr1.cs_ins - self.ctr1.cs_del,
            self.ctr1.cs_ins, self.ctr1.cs_ins_avg.avg(), self.ctr1.cs_del, self.ctr1.cs_del_avg.avg(),
            self.ctr2.nhits, self.ctr2.cs_good, self.ctr2.cs_good - self.ctr2.cs_ins - self.ctr2.cs_del,
            self.ctr2.cs_ins, self.ctr2.cs_ins_avg.avg(), self.ctr2.cs_del, self.ctr2.cs_del_avg.avg()))
        
    def writeMatrices(self):
        prefix = os.path.splitext(self.filename)[0]
        with open(prefix + ".mat1.csv", "w") as out:
            self.ctr1.writeMatrix(out)
        with open(prefix + ".mat2.csv", "w") as out:
            self.ctr2.writeMatrix(out)

class Main(object):
    infiles = []
    outfile = "results"
    perc = False

    nfiles = 0
    results = None

    def __init__(self):
        self.results = [ [], [], [], [] ]

    def parseArgs(self, args):
        if "-h" in args or "--help" in args:
            return False
        prev = ""
        for a in args:
            if prev == "-o":
                self.outfile = a
                prev = ""
            elif a in ["-o"]:
                prev = a
            elif a == "-p":
                self.perc = True
            else:
                self.infiles.append(a)
        self.nfiles = len(self.infiles)
        return (self.nfiles > 0) # I don't like that 0 == False

    def usage(self):
        sys.stdout.write("""scan.py - Scan sequences from fastq files and compute counts of sgRNA occurrences.

Usage: scan.py [options] fastq...

Where options are:

  -o O | Use O as base name for output files (default: '{}').
  -p   | Output percentages where appropriate instead of counts.

Fastq files can be compressed with gzip. 

The program generates three sets of tables, each one composed of several sheets. The first sheet
is a description of the contents of the other sheets. Files are named according to the following pattern:

  O.tableN.sheetM.txt

where O is set with the -o option, N is the table number, and M is the sheet number. All files are 
tab-delimited.

In addition, the program also generates one matrix for each input file, containing the base frequencies for
all positions in sgRNA1. These files are named as follows:

  O.matN.txt

where N is the index of the input file.

""".format(self.outfile))

    def run(self):
        for arg in self.infiles:
            sys.stderr.write("{}".format(arg))
            S = Scanner(arg)
            self.results[0].append(arg)
            S.run3(self)
            self.results[1].append(S.nreads)
            sys.stderr.write("\t{}".format(S.nreads))
            for k in ['h1' ,'nh1', 'h2', 'nh2']:
                sys.stderr.write("\t{:.2f}".format(S.averagers[k].avg()))
            sys.stderr.write("\n")

        self.writeMatrices()
        self.writeResults()
        self.writeResultsInsDel()
        self.writeResultsCons()

    def writeMatrices(self):
        for i in range(self.nfiles):
            with open("{}.mat{}.txt".format(self.outfile, i+1), "w") as out:
                c = self.results[2][i]
                c.writeMatrix(out)

    def writeResults(self):
        # Table 1 - hits, combinations of substitutions, insertions, deletions.
        with open(self.outfile + ".table1.README.txt", "w") as out:
            out.write("""Sheets:
            
G1\tGuide 1
G2.1\tGuide 2, wildtype (no mutations)
G2.2\tGuide 2, single mutation
G2.3\tGuide 2, unobserved haplotype (ignore this!)
G2.4\tGuide 2, double mutation

Legend:

Reads\tTotal reads in fastq file
Hits\tReads having proper primer and target sequence
M\tPerfect match 
S\tSubstitutions present
I\tInsertions present only
D\tDeletions present only
ID\tInsertions and deletions of equal length
ID*\tInsertions and deletions, net increase in target length
ID**\tInsertions and deletions, net decrease in target length
""")
        with open(self.outfile + ".table1.sheet1.txt", "w") as out:
            out.write("Filename\tReads\tG1_Hits\t" + "\t".join([ "G1_"+m for m in COLNAMES ]) + "\n")
            for i in range(self.nfiles):
                c = self.results[2][i]
                out.write("{}\t{}\t{}".format(self.results[0][i], self.results[1][i], c.nhits))
                for x in COUNTERS:
                    if self.perc:
                        out.write("\t{}".format(1.0*c.counters[x]/c.nhits))
                    else:
                        out.write("\t" + str(c.counters[x]))
                out.write("\n")
            out.write("\n")

        for idx in range(4):
            with open("{}.table1.sheet{}.txt".format(self.outfile, idx+1), "w") as out:
                out.write("Filename\tReads\tG2.{}_Hits\t".format(idx+1) + "\t".join([ "G2.{}_{}".format(idx+1, m) for m in COLNAMES ]) + "\n")
                for i in range(self.nfiles):
                    c = self.results[3][i]
                    out.write("{}\t{}\t{}".format(self.results[0][i], self.results[1][i], c.nhits[idx]))
                    for x in COUNTERS:
                        if self.perc:
                            out.write("\t{}".format(1.0*c.counters[idx][x]/c.nhits[idx]))
                        else:
                            out.write("\t" + str(c.counters[idx][x]))
                    out.write("\n")
                out.write("\n")

    def writeResultsInsDel(self):
        # Table 2 - Insertions and deletions stats

        with open(self.outfile + ".table2.README.txt", "w") as out:
            out.write("""Sheets:

G1\tGuide 1
G2.1\tGuide 2, wildtype (no mutations)
G2.2\tGuide 2, single mutation
G2.3\tGuide 2, unobserved haplotype (ignore this!)
G2.4\tGuide 2, double mutation

Legend:

Reads\tTotal reads in fastq file
Hits\tReads having proper primer and target sequence
Longer\tNumber of reads containing longer than expected target
LongerAvg\tAverage net increase
I1\tReads with 1nt increase
I2\tReads with 2nt increase
I3\tReads with 3nt increase
I4\tReads with 4nt increase
I5\tReads with 5nt increase
I+5\tReads with increase over 5nt
Shorter\tNumber of reads containing shorter than expected target
ShorterAvg\tAverage net decrease
D1\tReads with 1nt decrease
D2\tReads with 2nt decrease
D3\tReads with 3nt decrease
D4\tReads with 4nt decrease
D5\tReads with 5nt decrease
D+5\tReads with decrease over 5nt
""")

        with open(self.outfile + ".table2.sheet1.txt", "w") as out:
            out.write("Filename\tReads\tG1_Hits\t" + "\t".join([ "G1_Longer",  "G1_LongerAvg", "G1_I1", "G1_I2", "G1_I3", "G1_I4", "G1_I5", "G1_I+5",
                                                                 "G1_Shorter", "G1_ShorterAvg", "G1_D1", "G1_D2", "G1_D3", "G1_D4", "G1_D5", "G1_D+5" ]) + "\n")
            for i in range(self.nfiles):
                c = self.results[2][i]
                totlonger = c.counters[2] + c.counters[3] + c.counters[14] + c.counters[15]
                totshorter = c.counters[4] + c.counters[5] + c.counters[22] + c.counters[23]
                if self.perc:
                    out.write("{}\t{}\t{}\t{}\t{}".format(self.results[0][i], self.results[1][i], c.nhits, 1.0*totlonger/c.nhits, c.insavg.avg()))
                else:
                    out.write("{}\t{}\t{}\t{}\t{}".format(self.results[0][i], self.results[1][i], c.nhits, totlonger, c.insavg.avg()))
                for x in [1, 2, 3, 4, 5, 6]:
                    if self.perc:
                        out.write("\t{}".format(1.0*c.inslen[x]/c.nhits))
                    else:
                        out.write("\t{}".format(c.inslen[x]))
                if self.perc:
                    out.write("\t{}\t{}".format(1.0*totshorter/c.nhits, c.delavg.avg()))
                else:
                    out.write("\t{}\t{}".format(totshorter, c.delavg.avg()))
                for x in [1, 2, 3, 4, 5, 6]:
                    if self.perc:
                        out.write("\t{}".format(1.0*c.dellen[x]/c.nhits))
                    else:
                        out.write("\t{}".format(c.dellen[x]))
                out.write("\n")
            out.write("\n")

        for idx in range(4):
            with open("{}.table2.sheet{}.txt".format(self.outfile, idx+1), "w") as out:
                out.write("Filename\tReads\tG2.{idx}_Hits\tG2.{idx}_Longer\tG2.{idx}_LongerAvg\tG2.{idx}_I1\tG2.{idx}_I2\tG2.{idx}_I3\tG2.{idx}_I4\tG2.{idx}_I5\tG2.{idx}_I+5\tG2.{idx}_Shorter\tG2.{idx}_ShorterAvg\tG2.{idx}_D1\tG2.{idx}_D2\tG2.{idx}_D3\tG2.{idx}_D4\tG2.{idx}_D5\tG2.{idx}_D+5\n".format(**{'idx': idx+1}))
                for i in range(self.nfiles):
                    c = self.results[3][i]
                    totlonger = c.counters[idx][2] + c.counters[idx][3] + c.counters[idx][14] + c.counters[idx][15]
                    totshorter = c.counters[idx][4] + c.counters[idx][5] + c.counters[idx][22] + c.counters[idx][23]
                    if self.perc:
                        out.write("{}\t{}\t{}\t{}\t{}".format(self.results[0][i], self.results[1][i], c.nhits[idx], 1.0*totlonger/c.nhits[idx], c.insavg[idx].avg()))
                    else:
                        out.write("{}\t{}\t{}\t{}\t{}".format(self.results[0][i], self.results[1][i], c.nhits[idx], totlonger, c.insavg[idx].avg()))
                    for x in [1, 2, 3, 4, 5, 6]:
                        if self.perc:
                            out.write("\t{}".format(1.0*c.inslen[idx][x]/c.nhits[idx]))
                        else:
                            out.write("\t{}".format(c.inslen[idx][x]))
                    if self.perc:
                        out.write("\t{}\t{}".format(1.0*totshorter/c.nhits[idx], c.delavg[idx].avg()))
                    else:
                        out.write("\t{}\t{}".format(totshorter, c.delavg[idx].avg()))
                    for x in [1, 2, 3, 4, 5, 6]:
                        if self.perc:
                            out.write("\t{}".format(1.0*c.dellen[idx][x]/c.nhits[idx]))
                        else:
                            out.write("\t{}".format(c.dellen[idx][x]))
                    out.write("\n")
                out.write("\n")

    def writeResultsCons(self):
        with open(self.outfile + ".table3.README.txt", "w") as out:
            out.write("""Legend:
Reads\tTotal number of reads in fastq file
NoHits1\tNumber of reads with primer1 not containing the sgRNA1 sequence
Cons1\tNumber of reads in NoHits1 containing conserved sequence 1
Expected1\tNumber of reads with conserved sequence 1 in the expected position
Ins1\tNumber of reads with insertions between primer and conserved sequence
AvgIns1\tAverage insertion length
Del1\tNumber of reads with deletions between primer and conserved sequence
AvgDel1\tAverage deletion length

The next 7 columns are the same, for reads containing the primer2 sequence (3' end of the amplicon).
""")
        with open(self.outfile + ".table3.sheet1.txt", "w") as out:
            out.write("Filename\tReads\tNoHits1\tCons1\tExpected1\tIns1\tAvgIns1\tDel1\tAvgDel1\tNoHits2\tCons\tExpected2\tIns2\tAvgIns2\tDel2\tAvgDel2\n")
            for i in range(self.nfiles):
                ctr1 = self.results[2][i]
                ctr2 = self.results[3][i]
                out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                    self.results[0][i], self.results[1][i],
                    ctr1.nhits, ctr1.cs_good, ctr1.cs_good - ctr1.cs_ins - ctr1.cs_del,
                    ctr1.cs_ins, ctr1.cs_ins_avg.avg(), ctr1.cs_del, ctr1.cs_del_avg.avg(),
                    sum(ctr2.nhits), ctr2.cs_good, ctr2.cs_good - ctr2.cs_ins - ctr2.cs_del,
                    ctr2.cs_ins, ctr2.cs_ins_avg.avg(), ctr2.cs_del, ctr2.cs_del_avg.avg()))


if __name__ == "__main__":
    args = sys.argv[1:]
    M = Main()
    if M.parseArgs(args):
        M.run()
    else:
        M.usage()
