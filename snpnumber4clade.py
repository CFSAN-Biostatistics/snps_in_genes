#!/usr/bin/env python
#find SNP shared within a clade associated with the Agona manuscript comparing outbreaks between 1998 and 2008
#input is specific to that study and results files

from Bio import SeqIO

inputDir = "./"
textInputFile = open(inputDir + "snpma_preserved_with_reference.fasta","r")
itemFile = open(inputDir + "2008group.txt","r")
textOutputFile = open(inputDir + "output_2008group.txt","w")

#location found
locations = []
seqHash = dict()
ids = itemFile.readline().split()
#print ids[0],ids[1],ids[3]
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

textOutputFile.write(" ".join(locations))
print " ".join(locations)
textInputFile.close()
itemFile.close()
textOutputFile.close()
