from fastcomp import compare
import json

class ReadCompressor:
    def __init__(self):
        self.ampliconCountDict = {}
        self.totalReads = 0
        self.compressorAllowedThreshold = 500
        self.compressorMergeThreshold = 50
        with open('config.json') as config_file: 
            config_data = json.load(config_file)
            if ( 'j3x_MinReadsForTemplate' in config_data ):
                self.compressorAllowedThreshold = config_data['j3x_MinReadsForTemplate']
            if ( 'j3x_MaxReadsForMerge' in config_data ):
                self.compressorMergeThreshold = config_data['j3x_MaxReadsForMerge']
            
        if ( self.compressorAllowedThreshold <= self.compressorMergeThreshold ):
            print("ERROR: compressorAllowedThreshold must be greater than self.compressorMergeThreshold")
        
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
    
    def getDataList(self):
        readTupleList = list(self.ampliconCountDict.items())
        onesList = [x for x in readTupleList if x[1][0] <= self.compressorMergeThreshold ]
        filteredTupleList = [x for x in readTupleList if x[1][0] >= self.compressorAllowedThreshold]
        filteredTupleList.sort(key=lambda tup: -tup[1][0])
        matchedOnes = 0
        matchedMoreThanOne = 0
        totalOnes = 0
        for ones in onesList:
            matchedTupleList = [] # Stores the list of tuples that are 1 or 2 distance from the ones
            onesOccurence = ones[1][0]
            totalOnes += onesOccurence
            for filteredTuple in filteredTupleList:
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
                
        return filteredTupleList, totalOnes, matchedOnes, matchedMoreThanOne    