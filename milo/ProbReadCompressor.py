"""
ProbReadCompressor.py
This contains the class that handles the probabilistic read compressor.
ProbReadCompressor is executed in a few steps.
1. Input: Paired Reads with 
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

class ProbReadCompressor:
    def __init__(self, referenceCount, configFile = 'config.json'):
        self.ampliconCountDict = {}
        self.totalReads = 0
        self.referenceCount = referenceCount

    def putPairedRead(self, data):
        self.totalReads += 1
        sequenceInfo = data.sequenceInfo
        seq = data.sequenceData
        quality = data.qualityData
        print(seq)
        print(quality)
        # if ( seq not in self.ampliconCountDict ):
        #     # [Total Count, count from merges, sequence info, quality hash]
        #     self.ampliconCountDict[seq] = [0, 0, sequenceInfo, quality]
        # self.ampliconCountDict[seq][0] += 1
        