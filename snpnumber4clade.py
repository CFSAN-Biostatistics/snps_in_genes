#!/usr/bin/env python
#find SNP shared within a clade

from Bio import SeqIO
from optparse import OptionParser
import sys
import os.path
import csv


### Setup options
usage="""Takes a fasta file and filters reads based on a provided minimum length.
usage: %prog -i FILE [options]"""


parser = OptionParser(usage=usage)


parser.add_option("-i", "--in-file", metavar="FILE", dest="in_file", default=None,
				help="Specify the input SNP matrix FILE")


parser.add_option("-g", "--group", metavar="FILE", dest="group_file", default=None,
				help="Specify the file containing the sample names for a given group. THEY SHOULD BE SPACE DELIMITED - NOT NEW LINE.")


parser.add_option("-o", "--output-file", metavar="FILE", dest="out_file", default=None,
				help="Specify the output file. Default is prefix of the input file followed by '.DiagnosticSNPs.txt'")

(options, args) = parser.parse_args()


### Input file options
if options.in_file:
	in_file = os.path.abspath(options.in_file)
else:	
	print "You must specify an in_file. Use '-h' for help."
	sys.exit()

(in_filePath, in_fileWholeName) = os.path.split(in_file)
(in_fileBase, in_fileExt) = os.path.splitext(in_fileWholeName)


###	 Group file options
if options.group_file:
	group_file = os.path.abspath(options.group_file)
else:
	print "You must specify an group_file. Use '-h' for help."
	sys.exit()

(group_filePath, group_fileWholeName) = os.path.split(group_file)
(group_fileBase, group_fileExt) = os.path.splitext(group_fileWholeName)

###	 Output options
if options.out_file:
	out_file = os.path.abspath(options.out_file)
else:
	out_file = os.path.join(group_filePath, group_fileBase + "DiagnosticSNPs.txt")
	out_file = os.path.abspath(out_file)	

(out_filePath, out_fileWholeName) = os.path.split(out_file)
(out_fileBase, out_fileExt) = os.path.splitext(out_fileWholeName)


textInputFile = open(in_file,"r")
itemFile = open(group_file,"r")
textOutputFile = open(out_file,"w")

#location found
locations = []
seqHash = dict()
ids = itemFile.readline().split()
seqRecords = list(SeqIO.parse(textInputFile,"fasta"))
for seqRecord in seqRecords:
    seqString = str(seqRecord.seq)
    id = str(seqRecord.id)
    seqHash[id] = seqString

i = 0
itemId = ids[0]
itemStr = seqHash[itemId]
maxLen = len(itemStr)

SeqCount = len(seqRecords)
while i < maxLen:
    matchFlag = True
    for seqRecord in seqRecords:
        seqString = str(seqRecord.seq)
        curId = str(seqRecord.id)
        if (curId not in ids and (seqString[i] == itemStr[i])) \
        or (curId in ids and ((seqString[i] != itemStr[i]) or (seqString[i] == '-'))):
            matchFlag = False
            break

    if matchFlag:
        locations.append(str(i))

    i +=1

for i in locations:
	textOutputFile.write(group_fileBase + "\t" + "position" + str(int(i) + 1) + "\n")
textInputFile.close()
itemFile.close()
textOutputFile.close()
