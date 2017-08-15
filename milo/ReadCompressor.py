from genomicsUtils import simpleDistance

class ReadCompressor:
    def __init__(self):
        self.ampliconCountDict = {}
        self.totalReads = 0
        
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
        onesList = [x for x in readTupleList if x[1][0] <= 1 ]
        filteredTupleList = [x for x in readTupleList if x[1][0] > 2]
        filteredTupleList.sort(key=lambda tup: -tup[1][0])
        
        for ones in onesList:
            for tuples in filteredTupleList:
                dist = simpleDistance(ones[0], tuples[0])
                if ( dist < 2 ):
                    # print("Match\n{0}\n{1}".format(ones[0], tuples[0]))
                    tuples[1][0] += ones[1][0] # Add to total count
                    tuples[1][1] += ones[1][0] # Add to close match count
                    break
        return filteredTupleList