from j4xUtils import *

class MutationFinder:
        
    def __init__(self):
        pass
        
    def reinit(self):
        # Format of ampliconMutationHashDict
        # {
        #    'MUTATION_HASH': OCCURENCE_COUNT,
        #    'MUTATION_HASH': OCCURENCE_COUNT,
        #    ...
        # }
        self.ampliconMutationHashDict = {}
        
        self.referenceCount = 0
        self.ampliconRefs = []
        with open('references/Manifest.csv') as references:
            lineno = -1
            
            for line in references:
                if lineno >= 0:
                    self.referenceCount += 1
                    csvCells = line.split(',')
                    sequence = csvCells[3]
                    coordinates = int(csvCells[4])
                    
                    # ampliconRefs format
                    # [
                    #    [ SEQUENCE, COORDINATES, COUNT ]
                    # ]
                    self.ampliconRefs.append([sequence, coordinates, 0]);
                lineno += 1
    
    def getReferenceAmpliconArray(self):
        return self.ampliconRefs
        
    def getReferenceAmplicon(self, ampliconID):
        if ( ampliconID == 0 ):
            return None
        return self.ampliconRefs[ampliconID-1]
        
    def getMutationHashDict(self):
        return self.ampliconMutationHashDict
    
    # Puts the identified mutation into the hash
    def putMutationHash(self, ampliconID, mutationHash, referenceCoordinate):
        if ( ampliconID == 0 ):
            return
        key = "{0} {1}".format(ampliconID, mutationHash)
        referenceAmplicon = self.getReferenceAmplicon(ampliconID)
        referenceAmplicon[2] += 1 # Increment the count
        
        if ( len(mutationHash) == 0 ):
            return
        if ( key not in self.ampliconMutationHashDict ):
            self.ampliconMutationHashDict[key] = 0
        self.ampliconMutationHashDict[key] += 1
    
    def extractHighestOccuringMutations(self, minOccurences):
        mutationTupleList = list(self.getMutationHashDict().items())
        filteredTupleList = [x for x in mutationTupleList if x[1] >= minOccurences]
        filteredTupleList.sort(key=lambda tup: -tup[1])
        return filteredTupleList
        
    def identifyMutations(self, data):
        ampliconID = int(data[0][3:data[0].index(',')])
        if ( ampliconID == 0 ):
            return None, None, None
            
        sequenceData = data[1][:-1]
        referenceAmplicon = self.getReferenceAmplicon(ampliconID)
        
        referenceSequence = referenceAmplicon[0]
        referenceCoordinate = referenceAmplicon[1]
        mutationHash = mutationIDAsHash(referenceSequence, sequenceData)
        
        return ampliconID, mutationHash, referenceCoordinate
