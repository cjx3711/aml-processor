"""
This file contains the lookup table and helpful functions when
converting from ampID to geneID and stuff like that.
"""

from DataTypes import *

class GeneMap:
    def __init__(self):
        # Stores the gene manifest
        self.geneManifest = []
        # Maps amplicons to genes
        self.ampToGene = []
        
        with open('references/GeneManifest.csv') as geneManifestFile:
            next(geneManifestFile)
            
            maxAmpID = 0
            for line in geneManifestFile:
                line = line.split(',')
                geneID = int(line[0])
                count = int(line[1])
                ampIDs = [ int(x) for x in line[2].split(' ') ]
                name = line[3]
                sequence = line[4]
                firstCoord = int(line[5])
                lastCoord = int(line[6])
                startCoords = [ int(x) for x in line[7].split(' ') ]
                endCoords = [ int(x) for x in line[8].split(' ') ]
                
                for ampID in ampIDs:
                    maxAmpID = max(ampID, maxAmpID)
                
                singleGeneRef = GeneReference(geneID, count, ampIDs, name, sequence, firstCoord, lastCoord, startCoords, endCoords)
                self.geneManifest.append(singleGeneRef)
                
            self.ampToGene = [0] * (maxAmpID + 1)
            for geneRef in self.geneManifest:
                for ampID in geneRef.ampIDs:
                    ampID = int(ampID)
                    self.ampToGene[ampID] = geneRef.geneID
        
    def fromSameGene(ampID1, ampID2):
        ampID1, ampID2 = int(ampID1), int(ampID2)
        return self.ampToGene[ampID1] == self.ampToGene[ampID2]
    
    def nextToEachOther(ampID1, ampID2):
        pass
    
    def getCoordinates(ampID1, ampID2):
        pass
    
    def getSequence(ampID1, ampID2):
        pass

test = GeneMap()
                