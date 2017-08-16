from fastcomp import compare
import json

class ReadCompressor:
    def __init__(self):
        self.ampliconCountDict = {}
        self.totalReads = 0
        self.j3x_minVAFForTemplate = 0.05
        self.j3x_maxReadsForMerge = 50
        with open('config.json') as config_file: 
            config_data = json.load(config_file)
            if ( 'j3x_minVAFForTemplate' in config_data ):
                self.j3x_minVAFForTemplate = config_data['j3x_minVAFForTemplate']
            if ( 'j3x_maxReadsForMerge' in config_data ):
                self.j3x_maxReadsForMerge = config_data['j3x_maxReadsForMerge']
            if ( 'j3x_readDeletorThreshold' in config_data ):
                self.j3x_readDeletorThreshold = config_data['j3x_readDeletorThreshold']
            
        #if ( self.j3x_minVAFForTemplate <= self.j3x_maxReadsForMerge ):
        #    print("ERROR: j3x_minVAFForTemplate must be greater than j3x_maxReadsForMerge")
        
    def putPairedRead(self, data):
        self.totalReads += 1
        iddata = data[0]
        sequence = data[1]
        quality = data[2]
        key = sequence
        if ( key not in self.ampliconCountDict ):
            # [Total Count, Count from close match, ID data, Quality Hash]
            self.ampliconCountDict[key] = [0, 0, iddata, quality]
        self.ampliconCountDict[key][0] += 1
    
    def getRawDataList(self):
        readTupleList = list(self.ampliconCountDict.items())
        return readTupleList

    def getDataList(self, ampliconCounts):
        readTupleList = list(self.ampliconCountDict.items())
        onesList = []
        leftoverList = []
        templateTupleList = []
        for x in readTupleList:
            ampID = int(x[1][2].split(',')[0][3:])
            seqReadCount = x[1][0]
            if seqReadCount > self.j3x_maxReadsForMerge: # If number of reads is too high for merging into template, check if
                if (seqReadCount / ampliconCounts[ampID - 1]) >= self.j3x_minVAFForTemplate: # It qualifies for a template by having a high VA
                    templateTupleList.append(x)
                elif seqReadCount > self.j3x_readDeletorThreshold: # Otherwise, check to make sure the read count isn't high before discarding sequence
                    leftoverList.append(x)
            else:  # Otherwise, add it to the to-be-merged list
                onesList.append(x)
        templateTupleList.sort(key=lambda tup: -tup[1][0])
        matchedOnes = 0
        matchedMoreThanOne = 0
        totalOnes = 0
        for ones in onesList:
            matchedTupleList = [] # Stores the list of tuples that are 1 or 2 distance from the ones
            onesOccurence = ones[1][0]
            totalOnes += onesOccurence
            for filteredTuple in templateTupleList:
                dist = compare(ones[0], filteredTuple[0])
                if ( dist <= 2 and dist != -1 ):
                    # filteredTuple[1][0] += onesOccurence # Add to total count
                    # filteredTuple[1][1] += onesOccurence # Add to close match count
                    # matchedOnes += onesOccurence
                    matchedTupleList.append(filteredTuple)
                    break
            matchCount = len(matchedTupleList)
            if matchCount > 0:
                if matchCount > 1:
                    matchedMoreThanOne += onesOccurence
                matchedOnes += onesOccurence
                splitValue = onesOccurence / matchCount
                for matchedTuple in matchedTupleList:
                    filteredTuple[1][0] += splitValue
                    filteredTuple[1][1] += splitValue
            elif onesOccurence >= self.j3x_readDeletorThreshold:
                leftoverList.append(ones)
        # Merge the two lists
        filteredTupleList = templateTupleList + leftoverList
        return filteredTupleList, totalOnes, matchedOnes, matchedMoreThanOne    