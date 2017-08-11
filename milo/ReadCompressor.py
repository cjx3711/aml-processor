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
            self.ampliconCountDict[key] = [0, iddata, quality]
        self.ampliconCountDict[key][0] += 1
    
    def getDataList(self):
        readTupleList = list(self.ampliconCountDict.items())
        readTupleList.sort(key=lambda tup: -tup[1][0])
        return readTupleList