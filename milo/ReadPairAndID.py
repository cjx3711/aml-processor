from ReadPairer import *
from AmpliconMatcherHashSweep import *

"""
This does what the old ReadPairer used to do. It pairs unpaired reads and
assigns them to their respective IDs.

It makes use of ReadPairer to pair the reads and AmpliconMatcherHashSweep to ID the reads.

"""
class ReadPairAndID:
    def __init__(self, configFile = 'config.json', referenceFile = 'references/Manifest.csv'):
        self.readPairer = ReadPairer(configFile)
        self.ampliconMatcher = AmpliconMatcherHashSweep(referenceFile)
        
    def getReferenceCount(self):
        return self.ampliconMatcher.getReferenceCount()
    
    def alignAndMerge(self, left, right):
        bases, quality, collisions, score = self.readPairer.mergeUnpaired(left[1].rstrip(), reverseComplement(right[1].rstrip()), left[3].rstrip(), right[3].rstrip()[::-1])
        
        failedToPair = 1 if collisions == '?' else 0
        
        # Checks which amplicon a read belongs to, and whether it is a translocation
        ampID, ampIDTrans, matchType = self.ampliconMatcher.findAmplicon(bases)
        
        
        # Joins amplicon number, collision number, and coordinate as new ID, and appends bases and quality
        if ampIDTrans != None:
            IDPart = 'TL:{0}/{1}'.format(ampID, ampIDTrans)
        else:
            IDPart = 'ID:{0}'.format(ampID)
            
        readData = ", ".join((IDPart, 'C:'+collisions))
        return AlignedAndMerged(failedToPair, matchType, readData, bases, quality)