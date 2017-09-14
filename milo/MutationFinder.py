from j4xUtils import *
from genomicsUtils import reverseComplement
from TranslocatedBlockMatcher import *
import csv
from pprint import pprint
class MutationFinder:
    def __init__(self, references = 'references/Manifest.csv'):
        self.references = references
        self.translocatedBlockMatcher = TranslocatedBlockMatcher()
        
    def reinit(self):
        # Format of ampliconMutationHashList
        # [
        #    ('MUTATION_HASH', OCCURENCE_COUNT),
        #    ('MUTATION_HASH', OCCURENCE_COUNT),
        #    ...
        # ]
        self.ampliconMutationHashList = []
        self.ampliconTranslocationList = []
        
        self.referenceCount = 0
        self.ampliconRefs = []
        with open(self.references) as refFile:
            # ampliconRefs format
            # [
            #    [ SEQUENCE, COORDINATES, READ_COUNT, MUTATION_COUNT, TRANSLOCATION_COUNT ]
            # ]
            self.ampliconRefs = [[reverseComplement(line[3]) if line[2] == "-" else line[3]] + 
                                [int(line[4])] + [0, 0, 0] for line in list(csv.reader(refFile, delimiter=','))[1:]]
            self.referenceCount = len(self.ampliconRefs)

    def getRefAmpliconArray(self):
        return self.ampliconRefs
        
    def getRefAmplicon(self, ampliconID):
        if ( ampliconID == 0 ):
            return None
        return self.ampliconRefs[ampliconID-1]
    
    # Puts the identified mutation into the hash
    def putMutationHash(self, ampliconID, mutationHash, referenceCoordinate, readCount):
        if ( ampliconID == 0 ):
            return
        key = "{0} {1}".format(ampliconID, mutationHash)
        
        # Increment the mutation count
        self.getRefAmplicon(ampliconID)[2] += readCount
        
        if ( len(mutationHash) == 0 ):
            return
            
        self.getRefAmplicon(ampliconID)[3] += readCount
        self.ampliconMutationHashList.append((key, readCount))
    
    def putTranslocationHash(self, ampID1, ampID2, matchingBlocks, readCount):
        if ( ampID1 == 0 or ampID2 == 0 ):
            return
        
        # Add to total count
        # self.getRefAmplicon(ampID1)[2] += readCount
        # self.getRefAmplicon(ampID2)[2] += readCount
        
        if ( len(matchingBlocks) == 0 ):
            return
        
        # Add to translocation count
        self.getRefAmplicon(ampID1)[4] += readCount
        self.getRefAmplicon(ampID2)[4] += readCount
        matchingBlockString = str(matchingBlocks)
        matchingBlockString = matchingBlockString.replace('Match', '').replace("'R1'", str(ampID1)).replace("'R2'", str(ampID2))
        self.ampliconTranslocationList.append(('{0} {1} {2}'.format(ampID1, ampID2, matchingBlockString), readCount))
        
    def filterAndSortMutants(self, minReads, minVAF, sampleAmpStats, translocations):
        if translocations:
            ampliconList = self.ampliconTranslocationList
        else:
            ampliconList = self.ampliconMutationHashList
        filteredTupleList = []
        # Filters out templates below read and VAF thresholds
        for x in ampliconList:
            ampID = int(x[0][:x[0].find(' ')])
            templateReads = x[1]
            totalReads = self.getRefAmplicon(ampID)[2]
            #totalReads = sampleAmpStats[ampID].totalReads
            if templateReads >= minReads and templateReads / totalReads >= minVAF:
                filteredTupleList.append(x)
        # Sort mutants by read count
        filteredTupleList.sort(key=lambda tup: -tup[1])
        return filteredTupleList
        
    def identifyTranslocations(self, data):
        iddataParts = data[0].split(', ')
        amplicons = iddataParts[0][3:].split('/')
        ampID1 = int(amplicons[0])
        ampID2 = int(amplicons[1])
        if ampID1 == 0 or ampID2 == 0:
            return 'T', None, None, None, None
            
        readCount = int(float(iddataParts[2].strip()[2:]))
        sequenceData = data[1][:-1]
        
        refAmplicon1 = self.getRefAmplicon(ampID1)[0]
        refAmplicon2 = self.getRefAmplicon(ampID2)[0]
        
        matchingBlocks = self.translocatedBlockMatcher.findTranslocatedMatchingBlocks(sequenceData, refAmplicon1, refAmplicon2)
        
        return 'T', ampID1, ampID2, matchingBlocks, readCount
        
    def identifyMutations(self, data):
        # data[0] format:
        # ID:003, C:0, R:5, M:0
    
        # Handle translocations
        if ( data[0].startswith('TL:') ):
            return self.identifyTranslocations(data)
        
        iddataParts = data[0].split(', ')
        
        ampliconID = int(float(iddataParts[0][3:]))
        if ( ampliconID == 0 ):
            return 'M', None, None, None, None
        readCount = int(float(iddataParts[2].strip()[2:]))
        sequenceData = data[1][:-1]
        referenceAmplicon = self.getRefAmplicon(ampliconID)
        
        referenceSequence = referenceAmplicon[0]
        referenceCoordinate = referenceAmplicon[1]
        mutationHash = mutationIDAsHash(referenceSequence, sequenceData)
        
        return 'M', ampliconID, mutationHash, referenceCoordinate, readCount
