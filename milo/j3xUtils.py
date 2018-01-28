"""
All routines of the j3x file format.
"""

from genomicsUtils import *
from phredUtils import *

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
    qualityScores = "".join(simplifyQuality(phredToAccuDict[x]) for x in sequenceData[3][:-1])
    baseSeq = "".join("_" if x == "N" else x for x in sequenceData[1][:-1])
    return sequenceIdentifier + "\n" + baseSeq + "\n" + qualityScores + "\n" + "\n"
    
def extractAmpID(dataString):
    # Data string format:
    # ID:423, C:0, ...
    ampIDPart = dataString.split(',')[0].strip() # Format 'ID:423' or 'TL:123/456'
    ampID1 = int(ampIDPart[3:6])
    ampID2 = int(ampIDPart[7:10]) if ampIDPart.startswith('TL') else 0
    return ampID1, ampID2
