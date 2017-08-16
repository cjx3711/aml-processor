from fastcomp import compare
import json

class ReadCompressor:
    def __init__(self):
        self.ampliconCountDict = {}
        self.totalReads = 0
        self.j3x_minReadsForTemplate = 500
        self.j3x_maxReadsForMerge = 50
        with open('config.json') as config_file: 
            config_data = json.load(config_file)
            if ( 'j3x_minReadsForTemplate' in config_data ):
                self.j3x_minReadsForTemplate = config_data['j3x_minReadsForTemplate']
            if ( 'j3x_maxReadsForMerge' in config_data ):
                self.j3x_maxReadsForMerge = config_data['j3x_maxReadsForMerge']
            if ( 'j3x_readDeletorThreshold' in config_data ):
                self.j3x_readDeletorThreshold = config_data['j3x_readDeletorThreshold']
            
        if ( self.j3x_minReadsForTemplate <= self.j3x_maxReadsForMerge ):
            print("ERROR: j3x_minReadsForTemplate must be greater than j3x_maxReadsForMerge")
        
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
    def getDataList(self):
        readTupleList = list(self.ampliconCountDict.items())
        onesList = [x for x in readTupleList if x[1][0] <= self.j3x_maxReadsForMerge ]
        # Leftovers that should not be deleted. There may be some useful insights in there.
        leftoverList = [x for x in readTupleList if x[1][0] > self.j3x_maxReadsForMerge and x[1][0] < self.j3x_minReadsForTemplate and x[1][0] >= self.j3x_readDeletorThreshold ]
        templateTupleList = [x for x in readTupleList if x[1][0] >= self.j3x_minReadsForTemplate]
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