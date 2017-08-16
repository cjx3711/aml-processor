from genomicsUtils import simpleDistance
import json

class ReadCompressor:
    def __init__(self):
        self.ampliconCountDict = {}
        self.totalReads = 0
        
        with open('config.json') as config_file: 
            config_data = json.load(config_file)
            if ( 'compressorAllowedThreshold' in config_data ):
                self.compressorAllowedThreshold = config_data['compressorAllowedThreshold']
            if ( 'compressorMergeThreshold' in config_data ):
                self.compressorMergeThreshold = config_data['compressorMergeThreshold']
            
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
        totalOnes = 0
        for ones in onesList:
            onesOccurence = ones[1][0]
            totalOnes += onesOccurence
            for tuples in filteredTupleList:
                dist = simpleDistance(ones[0], tuples[0])
                if ( dist < 2 ):
                    # print("Match\n{0}\n{1}".format(ones[0], tuples[0]))
                    tuples[1][0] += onesOccurence # Add to total count
                    tuples[1][1] += onesOccurence # Add to close match count
                    matchedOnes += onesOccurence
                    break
                
        return filteredTupleList, totalOnes, matchedOnes        