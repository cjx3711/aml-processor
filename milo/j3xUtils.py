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