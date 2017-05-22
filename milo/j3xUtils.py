"""
All routines of the j3x file format.
"""

from genomicsUtils import *

qualityDict = getPhredQualityDict();
errorDict = getPhredErrorDict();

def simplifyQuality(p):
    """
    Returns the appropriate j3x quality score given probability p
    """
    if p >= 0.999:
        return "."
    elif p >= 0.99:
        return ","
    elif p >= 0.9:
        return "?"
    else:
        return "x"

def fastqToJ3X(sequenceData):
    """
    Converts a 4-element sequenceData, representing a 4-line FastQ sequence entry to its j3x equivalent
    """
    coordIndices = nthAndKthLetter(sequenceData[0], ":", 5, 7)
    sequenceIdentifier = sequenceData[0][coordIndices[0]: coordIndices[1] - 2]
    qualityScores = "".join(simplifyQuality(qualityDict[x]) for x in sequenceData[3][:-1])
    baseSeq = "".join("_" if x == "N" else x for x in sequenceData[1][:-1])
    return sequenceIdentifier + "\n" + baseSeq + "\n" + qualityScores + "\n" + "\n"
    

def getVal(hash, key):
    if key in hash:
        return hash[key]
    return 0
    
def quickDistance(sequence1, sequence2):
    hash1 = quickDistance(sequence1)
    hash2 = quickDistance(sequence2)
    dist = 0
    dist = dist + abs(len(hash1) - len(hash2))
    for i in range(min(len(hash1), len(hash2))):
        h1 = hash1[i]
        h2 = hash1[i]
        dist = abs(getVal(h1, 'A'), getVal(h2, 'A')) + abs(getVal(h1, 'T'), getVal(h2, 'T')) + abs(getVal(h1, 'C'), getVal(h2, 'C')) + abs(getVal(h1, 'G'), getVal(h2, 'G'))
    
    return dist
    
def quickHash(sequence):
    hash = []
    for i in range(len(a) // 4):
        partHash = {}
        seq = sequence[i * 4: (i+1) * 4]
        for c in seq:
            if c in partHash:
                partHash[c] = partHash[c] + 1
            else:
                partHash[c] = 1
        hash.append(partHash) 
    return hash
    
a = "ATCGCTACG"
b = "ATCGCTACG"
print(quickDistance(a,b))