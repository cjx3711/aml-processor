from j4xUtils import *

class MutationFinder:
        
    def __init__(self, groupsOf):
        self.groupsOf = groupsOf
        self.referenceCount = 0;
        self.ampliconRefs = []
        self.ampliconMutationMaps = {}
        with open('references/Manifest.csv') as references:
            lineno = -1
            
            for line in references:
                if lineno >= 0:
                    self.referenceCount += 1
                    sequence = line.split(',')[2]
                    
                    self.ampliconRefs.append(sequence);
                lineno += 1
        
    def getMutationMap(self):
        return self.ampliconMutationMaps
        
    def putMutationMap(self, ampliconID, mutationHash):
        key = "{0} {1}".format(ampliconID, mutationHash)
        if ( len(mutationHash) == 0 ):
            return
        if ( key not in self.ampliconMutationMaps ):
            self.ampliconMutationMaps[key] = 0
        self.ampliconMutationMaps[key] += 1
    
    def extractHighestOccuringMutations(self, minOccurences):
        mutationTupleList = list(self.getMutationMap().items())
        filteredTupleList = [x for x in mutationTupleList if x[1] >= minOccurences]
        filteredTupleList.sort(key=lambda tup: -tup[1])
        return filteredTupleList
        
    def identifyMutations(self, data):
        for i in range(self.groupsOf):
            ampliconID = int(data[4*i+0][3:data[4*i+0].index(',')])
            
            sequenceData = data[4*i+1][:-1]
            referenceSequence = self.ampliconRefs[ampliconID-1]
            mutationHash = mutationIDAsHash(referenceSequence, sequenceData)
            
            self.putMutationMap(ampliconID, mutationHash)
        
        return 0
