from fastcomp import compare
import json
from tqdm import tqdm
from j3xUtils import *

"""
DEPRECATED: Replaced with ProbReadCompressor
"""

class ReadCompressor:
    def __init__(self, referenceCount, configFile = 'config.json'):
        self.ampliconCountDict = {}
        self.totalReads = 0
        self.referenceCount = referenceCount
        self.j3x_minVAFForTemplate = 0.05
        self.j3x_maxReadsForMerge = 10
        self.j3x_readDeletorThreshold = 5
        with open(configFile) as config_file: 
            config_data = json.load(config_file)
            if ( 'j3x_minVAFForTemplate' in config_data ):
                self.j3x_minVAFForTemplate = config_data['j3x_minVAFForTemplate']
            if ( 'j3x_maxReadsForMerge' in config_data ):
                self.j3x_maxReadsForMerge = config_data['j3x_maxReadsForMerge']
            if ( 'j3x_readDeletorThreshold' in config_data ):
                self.j3x_readDeletorThreshold = config_data['j3x_readDeletorThreshold']

    def putPairedRead(self, data):
        self.totalReads += 1
        sequenceInfo = data.sequenceInfo
        seq = data.sequenceData
        quality = data.qualityData
        if ( seq not in self.ampliconCountDict ):
            # [Total Count, count from merges, sequence info, quality hash]
            self.ampliconCountDict[seq] = [0, 0, sequenceInfo, quality]
        self.ampliconCountDict[seq][0] += 1
    
    def prepareForCompression(self):
        # Init variables for compressing step
        self.ampliconCountList = [0] * (self.referenceCount + 1)
        self.discardCountList = [0] * (self.referenceCount + 1)
        self.tbmergedList = [] # 1. Sequences to be merged with template
        self.leftoverList = [] # 2. Sequences that cannot be merged, or cannot qualify as a template, but still have a high-ish read depth
        self.discardedList = [] # 3. Store the discarded stuff for possible further processing
        self.templateNestedList = [ [] for x in range(self.referenceCount + 1)] # 3. Double-nested list containing 572 indices for the amplicons, with each inner list holding all templates for that amplicon
        self.templateFlatList = [] # Same as above, but as one big flat list (all amplicons together) for easier output
        self.numOfTemplatesPerAmp = []
        
        self.numMergeAttempts = 0
        self.mergedCount = 0 # Only successful merges, both sure and unsure
        self.mergedUnsureCount = 0 # Count of merges where multiple candidate templates are equal distance away from sequence
        self.mergedD1Count = 0 # Count of merges with distance of 1
        self.mergedD2Count = 0 # and distance of 2
        self.failedMergeAndDiscarded = 0 # Count the number of items that were discarded after failing to merge 
        
        readTupleList = list(self.ampliconCountDict.items())
        # Computes the total read count for each amplicon to calculate VAF
        for seq in readTupleList:
            ampID1, ampID2 = extractAmpID(seq[1][2])
            seqReadCount = seq[1][0]
            self.ampliconCountList[ampID1] += seqReadCount
            self.ampliconCountList[ampID2] += seqReadCount if ampID2 != 0 else 0
                
        # Classify sequences into one of 3 groups, or discard them
        for seq in readTupleList:
            ampID1, ampID2 = extractAmpID(seq[1][2])
            seqReadCount = seq[1][0]
            if seqReadCount > self.j3x_maxReadsForMerge: # If number of reads is too high for merging into template, check if
                if seqReadCount / int(self.ampliconCountList[ampID1]) >= self.j3x_minVAFForTemplate: # It qualifies for a template by having a high VA
                    self.templateNestedList[ampID1].append(seq)
                    self.templateFlatList.append(seq)
                else:
                    if seqReadCount >= self.j3x_readDeletorThreshold: # Otherwise, check to make sure the read count isn't high before discarding sequence
                        self.leftoverList.append(seq)
                    else:
                        self.discardedList.append(seq)
                        self.discardCountList[ampID1] += 1 # Discard
            else:  # Otherwise, add it to the to-be-merged list
                self.tbmergedList.append(seq)

    def getStats(self):
        return ( self.numMergeAttempts,
                 self.mergedCount,
                 self.mergedUnsureCount,
                 self.mergedD1Count,
                 self.mergedD2Count,
                 self.discardCountList,
                 self.ampliconCountList,
                 self.failedMergeAndDiscarded,
                 len(self.templateFlatList),
                  # List of templates for each amplicon
                 [len(templatesPerAmp) for templatesPerAmp in self.templateNestedList]
                 )
                 
    def getDiscardedList(self):
        return self.discardedList;

    def compress(self):
        # Merging to-be-merged list with template list
        for seq in tqdm(self.tbmergedList):
            mergeCandidatesD1 = [] # List containing templates that each sequence might be merged with distance of 1,
            mergeCandidatesD2 = [] # and distance of 2
            seqReadCount = seq[1][0]
            ampID1, ampID2 = extractAmpID(seq[1][2])

            self.numMergeAttempts += seqReadCount

            for template in self.templateNestedList[ampID1]: # Get edit distance between sequence and every applicable template
                dist = compare(seq[0], template[0])
                if dist != -1: # If distance is not more than 2, put template in consideration for merge
                    if dist == 1:
                        mergeCandidatesD1.append(template)
                    else:
                        mergeCandidatesD2.append(template)

            numCandidates = len(mergeCandidatesD1) + len(mergeCandidatesD2)
            if numCandidates > 0:
                self.mergedCount += seqReadCount
                if numCandidates > 1:
                    self.mergedUnsureCount += seqReadCount
                if mergeCandidatesD1: # Prioritize templates that are a distance of 1, rather than 2, from the current sequence
                    splitValue = seqReadCount / len(mergeCandidatesD1) # Allocate read count equally among templates equally similar to sequence
                    self.mergedD1Count += seqReadCount
                    for template in mergeCandidatesD1:
                        template[1][0] += splitValue # Increase total read count
                        template[1][1] += splitValue # Increase read count of merges
                else:
                    splitValue = seqReadCount / len(mergeCandidatesD2)
                    self.mergedD2Count += seqReadCount
                    for template in mergeCandidatesD2:
                        template[1][0] += splitValue
                        template[1][1] += splitValue

            else:
                if seqReadCount >= self.j3x_readDeletorThreshold: # If we can't merge the sequence but it has a high read depth
                    self.leftoverList.append(seq)
                else:
                    self.discardedList.append(seq)
                    self.discardCountList[ampID1] += 1 # Discard
                    self.failedMergeAndDiscarded += 1
        # Combine the newly reinforced templates with the leftovers for inclusion in j3x
        j3xSeqs = self.templateFlatList + self.leftoverList
        return j3xSeqs