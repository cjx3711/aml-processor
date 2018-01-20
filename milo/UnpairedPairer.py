from collections import defaultdict
from genomicsUtils import grouper
from ReadPairer import *
from generalUtils import getFileList
from pprint import pprint
from GeneMap import *

"""
Opens a j3x file and attempts to pair the unpaired reads.

All existing reads will be output into another j3x file ending with _TILED.j3x.
The unpaired reads will then be tiled with adjacent amplicons if any.
"""

class UnpairedPairer:
    def __init__(self, configFile = 'config.json'):
        self.unpairedDict = defaultdict(list)
        self.mergedList = []
        # READ JSON WIH SCORE THRESHOLD HERE
        self.scoreThreshold = 10
        
        self.readPairer = ReadPairer(configFile)
        self.geneMap = GeneMap()

    def idSeqToVals(self, idSeq):
        """
        Returns a generator of ampID, numCollisions, numReads, numMerges.
        Note: Returns values as strings due to possibiliy of non-numeric in numCollisions.
        """
        return (''.join(list(filter(lambda x: str.isdigit(x) or x == '?', chunk))) for chunk in idSeq.split(", "))

    def populateUnpaired(self, j3xFilepath, tiledj3xFilepath):
        """
        Extracts all unpaired templates from compressed j3x, and deletes extracted lines from said j3x.
        Run this before pairUnpaired.
        """
        j3xInput = open(j3xFilepath)
        j3xOutput = open(tiledj3xFilepath, "w+", newline = "")
        
        for idSeq, bases, qualitySeq, _ in grouper(j3xInput, 4):

            ampIDString, numCollisions, numReads, numMerges = self.idSeqToVals(idSeq)
            if len(ampIDString) == 3 and ampIDString != "000" and numCollisions == "?": # Eliminate all with unknown ampIDs
                left, right = bases.split(" ")
                lQuality, rQuality = qualitySeq.split(" ")
                # print("{0} {1} {2}".format(ampIDString, numReads, bases))
                self.unpairedDict[int(ampIDString)].append((left.strip(), right.strip(), lQuality.strip(), rQuality.strip(), idSeq))
            else: # Output the rest of the j3x file to the outputFile
                j3xOutput.write(idSeq)
                j3xOutput.write(bases)
                j3xOutput.write(qualitySeq)
                j3xOutput.write(_)
        
        j3xInput.close()
        j3xOutput.close()
        
    def pairUnpaired(self):
        """
        Pairs every template with any amplicon in overlapping tiles above a score threshold, loads them into mergedList, and returns for writing.
        """
        ampIDs = sorted(list(self.unpairedDict.keys()))
        for ampID in ampIDs: # Don't iterate through unpairedDict directly otherwise idSeqs in mergedList might be wrong way round.
            stillUnmergedList = []
            for left, right, lQuality, rQuality, idSeq in self.unpairedDict[ampID]:
                # Find overlapping tiles on left and right of current amplicon. None type for either or both if no overlapping tile.
                leftTileID, rightTileID = self.geneMap.adjacentTiles(ampID)
                merged = False
                if leftTileID:
                    seqsAndScores = [self.pair(left, leftTileRight, lQuality, leftTileRQ) + (leftTileIDSeq,) for 
                                     leftTileLeft, leftTileRight, leftTileLQ, leftTileRQ, leftTileIDSeq in self.unpairedDict[leftTileID]]
                    if self.pickAndMerge(seqsAndScores, ampID, idSeq):
                        print("Merged left")
                        merged = True 
                    
                if rightTileID:
                    seqsAndScores = [self.pair(rightTileLeft, right, rightTileLQ, rQuality) + (rightTileIDSeq,) for
                                     rightTileLeft, rightTileRight, rightTileLQ, rightTileRQ, rightTileIDSeq in self.unpairedDict[rightTileID]]
                    if self.pickAndMerge(seqsAndScores, ampID, idSeq):
                        merged = True 
                
                if not merged:
                    stillUnmergedList.append((left, right, lQuality, rQuality, idSeq))
                        
            self.unpairedDict[ampID] = stillUnmergedList
            
        pprint(len(self.mergedList))
        
        self.mergedList = self.compressMergedList()
        pprint(len(self.mergedList))
        pprint(self.mergedList)
        
        return None
    
    def pair(self, newLeft, newRight, lQuality, rQuality):
        return self.readPairer.mergeUnpaired(newLeft, newRight, lQuality, rQuality, True)

    def pickAndMerge(self, seqsAndScores, ampIDNumeric, idSeq):
        merged = False
        mergedCount = 0 # How many did we manage to merge with?
        for mergedSeq, mergedQuality, collisions, mergeScore, idSeq2 in seqsAndScores:
            if mergeScore > self.scoreThreshold:
                mergedCount += 1
        for mergedSeq, mergedQuality, collisions, mergeScore, idSeq2 in seqsAndScores:
            if mergeScore > self.scoreThreshold:
                # print ("Merge score: {0}".format(mergeScore))
                self.addToMergedList(mergedSeq, mergedQuality, ampIDNumeric, idSeq, idSeq2, collisions, mergedCount)
                merged = True
        return merged

    def addToMergedList(self, mergedSeq, mergedQuality, ampIDNumeric, idSeq1, idSeq2, collisions, mergedCount):
        # THINK OF NEW ID SEQ THAT MAKES SENSE
        # I recommend starting with MT (merged tiles): ampID1, ampID2. Then reusing our translocation finder to find mutants.
        # Once again, reads dunno how to count
        ampIDString1, numCollisions1, numReads1, numMerges1 = self.idSeqToVals(idSeq1)
        ampIDString2, numCollisions2, numReads2, numMerges2 = self.idSeqToVals(idSeq2)
        numReads1 = int(int(numReads1) / mergedCount)
        newIDSeq = "MT:{0}/{1}, C:{2}, R:{3}/{4}".format(ampIDString1, ampIDString2, collisions, numReads1, numReads2)
        self.mergedList.append((newIDSeq, mergedSeq, mergedQuality))
    
    
    def compressMergedList(self): # Compresses the merged list based on the ID and sequence.
        mergedHash = {}
        for merged in self.mergedList:
            idSeq = merged[0]
            mergedSeq = merged[1]
            idSeqParts = idSeq.split(',')
            amps = idSeqParts[0].strip()
            counts = idSeqParts[2].strip()[2:].split('/')
            count1 = int(counts[0])
            count2 = int(counts[1])
            key = amps + mergedSeq
            if not key in mergedHash:
                mergedHash[key] = [merged, 0, 0]
            mergedHash[key][1] += count1
            mergedHash[key][2] += count2
            
        newMergedList = []
        for mergedMerged in list(mergedHash.values()):
            merged = mergedMerged[0]
            count1 = mergedMerged[1]
            count2 = mergedMerged[2]
            idSeq = merged[0]
            mergedSeq = merged[1]
            mergedQuality = merged[2]
            newIDSeq = idSeq[:idSeq.find('R:')+2] + str(count1) + '/' + str(count2)
            newMergedList.append((newIDSeq, mergedSeq, mergedQuality))
        return newMergedList

if __name__ ==  "__main__":
    directory = 'data/2-paired/'
    filenameArray = getFileList('files.json')
    for filenames in filenameArray:
        unpairedPairer = UnpairedPairer()
        unpairedPairer.populateUnpaired(directory + filenames.paired, directory + filenames.tiled)
    
        unpairedPairer.pairUnpaired()
    
    
    
    
    