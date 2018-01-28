"""
ProbReadCompressor.py
This contains the class that handles the probabilistic read compressor.
ProbReadCompressor is executed in a few steps.
1. Input: PairedRead tuple
2. Convert the reads into ProbBase format.
3. Collapse all the ProbBase reads into their counts.
4. Perform probabilistic compression on the reduced set.
   a. Similar to the other compression, but it takes into account the read
      probabilities when making decisions on what to merge.
"""

from fastcomp import compare
import json
from tqdm import tqdm
from j3xUtils import *
from phredUtils import *
from pprint import pprint

class ProbReadCompressor:
    def __init__(self, referenceCount, configFile = 'config.json'):
        # Format of probaseDict:
        # {
        #   probaseKey: [
        #       count, baseSeq, totalPhredScores
        #   ]
        # }
        self.probaseDict = {}
        self.totalReads = 0
        self.referenceCount = referenceCount

    def putPairedRead(self, pairedRead):
        # pprint(pairedRead)
        self.totalReads += 1
        if not pairedRead.pairSuccess:
            return
        baseSeq = pairedRead.baseSeq
        phredSeq = pairedRead.phredSeq
        probaseKey = getProbaseFromBaseAndPhred(baseSeq, phredSeq)
        phrednp = phredSeqToNpArray(phredSeq)
        # pprint(probaseKey)
        if ( probaseKey not in self.probaseDict ):
            self.probaseDict[probaseKey] = [1, baseSeq, phrednp]
        else:
            # Add the count and all the scores.
            # Call normalisePhredScores after all reads are in to divide.
            self.probaseDict[probaseKey][0] += 1
            self.probaseDict[probaseKey][2] += phrednp
        
        # if ( seq not in self.ampliconCountDict ):
        #     # [Total Count, count from merges, sequence info, quality hash]
        #     self.ampliconCountDict[seq] = [0, 0, sequenceInfo, quality]
        # self.ampliconCountDict[seq][0] += 1
        
    def normalisePhredScoreTotal(self):
        print("Compressed Length: {0}".format(len(self.probaseDict.keys())))
        # Go through the probaseDict and divide all the phredScores by the counts
        for key in self.probaseDict.keys():
            self.probaseDict[key][2] /= self.probaseDict[key][0]
    
    def compareScore(sequence1, sequence2):
        return 0
    
    def prepareForCompression(self):
        self.normalisePhredScoreTotal()
        # ID the amplicons and put into an amplicon list.
        
        # For each amplicon,
        #   For each compressedRead
        #       If it can merge with anything in template list
        #           Distrubute counts evenly
        #       else
        #           Add to template list
            
