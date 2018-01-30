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
from AmpliconMatcherHashSweep import *
from pprint import pprint

class ProbReadCompressor:
    def __init__(self, referenceCount, configFile = 'config.json', referenceFile = 'references/Manifest.csv'):
        # Format of probaseDict:
        # {
        #   probaseKey: [
        #       count, baseSeq, totalPhredScores
        #   ]
        # }
        self.probaseDict = {}
        self.totalReads = 0
        self.referenceCount = referenceCount
        
        self.ampliconMatcher = AmpliconMatcherHashSweep(referenceFile)

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
    
    """
    Returns a score of the match
    Returns 0 if the reads are not same length
    """
    def computeScore(compressedRead1, compressedRead2):
        baseSeq1 = compressedRead1[1]
        baseSeq2 = compressedRead2[1]
        phrednp1 = compressedRead1[2]
        phrednp2 = compressedRead2[2]
        if len(baseSeq1) != len(baseSeq2):
            return 0
        
        score = 0
        for bp1, bp2, pnp1, pnp2 in baseSeq1, baseSeq2, phrednp1, phrednp2:
            pass
        
        return score
    
    """
    Prepares the lists and hashes required for compression
    Sorts the compressedReads into their respective amplicon bins
    """
    def prepareForCompression(self):
        self.normalisePhredScoreTotal()
        
        # 1. Double-nested list containing 572 indices for the amplicons, with each inner list holding all templates for that amplicon
        self.compressedReadsByAmpID = [ [] for x in range(self.referenceCount + 1)]
        self.translocatedCompressedReads = []
        
        # ID the amplicons and put into an amplicon list.
        for key in self.probaseDict.keys():
            # compressedRead format [count, baseSeq, phrednp]
            compressedRead = self.probaseDict[key]
            ampID, ampIDTrans, matchType = self.ampliconMatcher.findAmplicon(compressedRead[1])
            ampID = int(ampID)
            if ampIDTrans != None:
                self.translocatedCompressedReads.append(compressedRead)
            else:
                self.compressedReadsByAmpID[ampID].append(compressedRead)
    
    """
    Match up reads that are close in edit distance to bring down filesize.
    """
    def compress(self):
        templates = [] # Used to store current templates
        # For each amplicon,
        for ampCompressedReads in self.compressedReadsByAmpID:
            templates.clear()
        #   For each compressedRead
            for compressedRead in ampCompressedReads:
                for templateRead in templates:
                    pass
        #       If it can merge with anything in template list
        #           Distrubute counts evenly
        #       else
        #           Add to template list
            
# compareScore("ATCG", "AAAT")
