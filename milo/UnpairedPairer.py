from collections import defaultdict
from genomicsUtils import grouper
from ReadPairer import *

"""
Opens a j3x file and attempts to pair the unpaired reads.

All existing reads will be output into another j3x file ending with _TILED.j3x.
The unpaired reads will then be tiled with adjacent amplicons if any.
"""

class UnpairedPairer:
    def __init__(self, configFile = 'config.json'):
        self.unpairedDict = defaultdict(list)
        self.mergedDict = {}
        # READ JSON WIH SCORE THRESHOLD HERE
        
        self.readPairer = ReadPairer(configFile)

    def idSeqToVals(self, idSeq):
        """
        Returns a generator of ampID, numCollisions, numReads, numMerges.
        Note: Returns values as strings due to possibiliy of non-numeric in numCollisions.
        """
        return (''.join(list(filter(str.isdigit, chunk))) for chunk in idSeq.split(", "))

    def populateUnpaired(self, j3xFilepath, tiledj3xFilepath):
        """
        Extracts all unpaired templates from compressed j3x, and deletes extracted lines from said j3x.
        Run this before pairUnpaired.
        """
        j3xInput = open(j3xFilepath)
        j3xOutput = open(tiledj3xFilepath, "w+", newline = "")
        
        for idSeq, bases, qualitySeq, _ in grouper(j3xInput, 4):

            ampIDString, numCollisions, numReads, numMerges = self.idSeqToVals(idSeq)
            if ampIDString != "000" and numCollisions == "?": # Eliminate all with unknown ampIDs
                left, right = bases.split(" ")
                lQuality, rQuality = qualitySeq.split(" ")
                self.unpairedDict[int(ampIDString)].append((left, right, lQuality, rQuality, idSeq))
            else: # Output the rest of the j3x file to the outputFile
                j3xOutput.write(idSeq)
                j3xOutput.write(bases)
                j3xOutput.write(qualitySeq)
                j3xOutput.write(_)
        
        j3xInput.close()
        j3xOutput.close()


    def pairUnpaired(self, refLength):
        """
        Pairs every template with any amplicon in overlapping tiles above a score threshold, loads them into mergedDict, and returns for writing.
        """
        for ampIDNumeric in range(1, refLength + 1): # Don't iterate through unpairedDict directly otherwise idSeqs in mergedDict might be wrong way round.
            for left, right, lQuality, rQuality, idSeq in unpairedDict[i]:

                    # Find overlapping tiles on left and right of current amplicon. None type for either or both if no overlapping tile.
                    leftTileID, rightTileID = self.overlappingTiles[ampIDNumeric]
                    if leftTileID:
                        seqsAndScores = [self.pair(left, leftTileRight, lQuality, leftTileRQ) + (leftTileIDSeq,) for 
                                         leftTileLeft, leftTileRight, leftTileLQ, leftTileRQ, leftTileIDSeq in self.unpairedDict[leftTileID]]
                        self.pickAndMerge(seqsAndScores, ampIDNumeric, idSeq)
                    if rightTileID:
                        seqsAndScores = [self.pair(rightTileLeft, right, rightTileLQ, rQuality) + (rightTileIDSeq,) for
                                         rightTileLeft, rightTileRight, rightTileLQ, rightTileRQ, rightTileIDSeq in self.unpairedDict[rightTileID]]
                        self.pickAndMerge(seqsAndScores, ampIDNumeric, idSeq)

        return mergedDict

    def pair(self, newLeft, newRight, lQuality, rQuality):
        """
        
        """
        return readPairer.mergeUnpaired(newLeft, newRight, lQuality, rQuality)

    def pickAndMerge(self, seqsAndScores, ampIDNumeric, idSeq):
        for mergedSeq, mergedQuality, collisions, mergeScore, idSeq2 in seqsAndScores:
            if mergeScore > self.scoreThreshold:
                self.updateMerge(mergedSeq, mergedQuality, ampIDNumeric, idSeq, idSeq2, collisions)

    def updateMerge(self, mergedSeq, mergedQuality, ampIDNumeric, idSeq, idSeq2, collisions):
        self.unpairedDict[ampIDNumeric].remove((left, right, idSeq))
        # THINK OF NEW ID SEQ THAT MAKES SENSE
        # I recommend starting with MT (merged tiles): ampID1, ampID2. Then reusing our translocation finder to find mutants.
        # Once again, reads dunno how to count
        # newIDSeq = "MT:{0}/{1}, C:{2}, ..."
        self.mergedDict.append((newIDSeq, mergedSeq, mergedQuality))
        
if __name__ ==  "__main__":
    unpairedPairer = UnpairedPairer()
    
    inputFile = 'data/2-paired/SMALLTEST_AD01_S1_L001_PAIRED.j3x'
    outputFile = 'data/2-paired/SMALLTEST_AD01_S1_L001_TILED.j3x'
    unpairedPairer.populateUnpaired(inputFile, outputFile)
    
    
    
    
    
    
    