from j4xUtils import *
from pprint import pprint

class GeneAmpMap:
    def __init__(self):
        self.ampRefCount = -1 # Don't count the first line
        # Count the number of references there are
        with open('references/Manifest.csv') as references:
            for line in references:
                self.ampRefCount += 1
                
        self.initGeneMap()
        
    
    def initGeneMap(self):
        self.geneRefs = []
        self.ampIDToGeneID = [ 0 for x in range(self.ampRefCount + 1) ]
        with open('references/GeneManifest.csv') as geneReferences:
            for line in geneReferences:
                splitParts = line.split(',')
                geneID = int(splitParts[0])
                ampIDString = splitParts[1]
                ampIDs = ampIDString.split(' ')
                geneName = splitParts[2]
                geneSequence = splitParts[3]
                startCoord = int(splitParts[4])
                endCoord = int(splitParts[5])
                
                ampIDList = []
                # Create ampIDToGeneID map
                for ampID in ampIDs:
                    ampID = int(ampID)
                    ampIDList.append(ampID)
                    self.ampIDToGeneID[ampID] = geneID
            
                # Create Gene map
                self.geneRefs.append((geneID, ampIDList, geneName, geneSequence, startCoord, endCoord))
            print(self.ampIDToGeneID)
            pprint(self.geneRefs)
    
    def getGeneIDFromAmpID(ampID):
        ampID = int(ampID)
        if ( ampID == 0 or ampID >= self.ampRefCount ):
            return 0
        return self.ampIDToGeneID[ampID]
    
gene = GeneAmpMap()