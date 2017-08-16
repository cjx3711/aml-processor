from fastcomp import compare
import json

class ReadCompressor:
    def __init__(self, referenceCount):
        self.ampliconCountDict = {}
        self.totalReads = 0
        self.referenceCount = referenceCount
        self.j3x_minVAFForTemplate = 0.05
        self.j3x_maxReadsForMerge = 10
        self.j3x_readDeletorThreshold = 5
        with open('config.json') as config_file: 
            config_data = json.load(config_file)
            if ( 'j3x_minVAFForTemplate' in config_data ):
                self.j3x_minVAFForTemplate = config_data['j3x_minVAFForTemplate']
            if ( 'j3x_maxReadsForMerge' in config_data ):
                self.j3x_maxReadsForMerge = config_data['j3x_maxReadsForMerge']
            if ( 'j3x_readDeletorThreshold' in config_data ):
                self.j3x_readDeletorThreshold = config_data['j3x_readDeletorThreshold']
        
        
    def putPairedRead(self, data):
        self.totalReads += 1
        seqInfoLine = data[0]
        seq = data[1]
        quality = data[2]
        key = seq
        if ( key not in self.ampliconCountDict ):
            # [Total Count, count from merges with distance 1, sequence data, quality hash]
            self.ampliconCountDict[key] = [0, 0, seqInfoLine, quality]
        self.ampliconCountDict[key][0] += 1
    
    def prepareForCompression(self):
        # Init variables for compressing step
        self.ampliconCountList = [0] * self.referenceCount
        self.discardCountList = [0] * self.referenceCount
        self.tbmergedList = [] # 1. Sequences to be merged with template
        self.leftoverList = [] # 2. Sequences that cannot be merged, or cannot qualify as a template, but still have a high-ish read depth
        self.templateList = [] # 3. Variant templates with high read depth
        
        self.numMergeAttempts = 0
        self.mergedCount = 0 # Only successful merges, both sure and unsure
        self.mergedUnsureCount = 0 # Count of merges where multiple candidate templates are equal distance away from sequence
        self.mergedD1Count = 0 # Count of merges with distance of 1
        self.mergedD2Count = 0 # and distance of 2
        
        readTupleList = list(self.ampliconCountDict.items())
        for seq in readTupleList:
            ampID = int(seq[1][2].split(',')[0].strip()[3:]) # Extracts the ampID from the info line
            seqReadCount = seq[1][0]
            self.ampliconCountList[ampID] += seqReadCount
            
            # Classify sequences into one of 3 groups below, or discard them
            if seqReadCount > self.j3x_maxReadsForMerge: # If number of reads is too high for merging into template, check if
                if (seqReadCount / self.ampliconCountList[ampID]) >= self.j3x_minVAFForTemplate: # It qualifies for a template by having a high VA
                    self.templateList.append(seq)
                else:
                    if seqReadCount > self.j3x_readDeletorThreshold: # Otherwise, check to make sure the read count isn't high before discarding sequence
                        self.leftoverList.append(seq)
                    else:
                         self.discardCountList[ampID] += 1 # Discard
            else:  # Otherwise, add it to the to-be-merged list
                self.tbmergedList.append(seq)
                
    def getStats(self):
        return self.numMergeAttempts, self.mergedCount, self.mergedUnsureCount, self.mergedD1Count, self.mergedD2Count, self.discardCountList, self.ampliconCountList
    
    def mergeIntoTemplates(self, seq):
        mergeCandidatesD1 = [] # List containing templates that each sequence might be merged with distance of 1,
        mergeCandidatesD2 = [] # and distance of 2

        for template in enumerate(self.templateList): # Get edit distance between sequence and every applicable template
            dist = compare(seq[0], template[1][0])
            if dist != -1: # If distance is not more than 2, put template in consideration for merge
                if dist == 1:
                    mergeCandidatesD1.append(template[0])
                else:
                    mergeCandidatesD2.append(template[0])
                    
        return mergeCandidatesD1, mergeCandidatesD2, seq
        
    def addMergeCadidates(self, mergeCandidatesD1, mergeCandidatesD2, seq):
        seqReadCount = seq[1][0]
        self.numMergeAttempts += seqReadCount
        numCandidates = len(mergeCandidatesD1) + len(mergeCandidatesD2)
        if numCandidates > 0:
            if numCandidates > 1:
                self.mergedUnsureCount += seqReadCount
            self.mergedCount += seqReadCount

            if mergeCandidatesD1: # Prioritize templates that are a distance of 1, rather than 2, from the current sequence
                splitValue = seqReadCount / len(mergeCandidatesD1) # Allocate read count equally among templates equally similar to sequence
                self.mergedD1Count += seqReadCount
                for templateInd in mergeCandidatesD1:
                    self.templateList[templateInd][1][0] += splitValue # Increase total read count
                    self.templateList[templateInd][1][1] += splitValue # Increase read count of merges
            else:
                splitValue = seqReadCount / len(mergeCandidatesD2)
                self.mergedD2Count += seqReadCount
                for templateInd in mergeCandidatesD2:
                    self.templateList[templateInd][1][0] += splitValue
                    self.templateList[templateInd][1][1] += splitValue

        else:
            if seqReadCount >= self.j3x_readDeletorThreshold: # If we can't merge the sequence but it has a high read depth
                self.leftoverList.append(seq)
            else:
                self.discardCountList[int(seq[1][2].split(',')[0][3:])] += 1 # Discard
        
        
    def getTbmerged(self):
        return self.tbmergedList
    
    def getMergedDataList(self):
        return self.templateList + self.leftoverList
        
    def getDataList(self):
        # Merging to-be-merged list with template list
        for seq in self.tbmergedList:
            mergeCandidatesD1 = [] # List containing templates that each sequence might be merged with distance of 1,
            mergeCandidatesD2 = [] # and distance of 2
            
            for template in self.templateList: # Get edit distance between sequence and every applicable template
                dist = compare(seq[0], template[0])
                if dist != -1: # If distance is not more than 2, put template in consideration for merge
                    if dist == 1:
                        mergeCandidatesD1.append(template)
                    else:
                        mergeCandidatesD2.append(template)
                        
            seqReadCount = seq[1][0]
            self.numMergeAttempts += seqReadCount
            numCandidates = len(mergeCandidatesD1) + len(mergeCandidatesD2)
            if numCandidates > 0:
                if numCandidates > 1:
                    self.mergedUnsureCount += seqReadCount
                self.mergedCount += seqReadCount

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
                    self.discardCountList[int(seq[1][2].split(',')[0][3:])] += 1 # Discard
        # Combine the newly reinforced templates with the leftovers for inckusion in j3x
        j3xSeqs = self.templateList + self.leftoverList
        return j3xSeqs