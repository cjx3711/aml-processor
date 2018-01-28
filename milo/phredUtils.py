from math import pow, log10
import numpy as np

phredSymbolList = ['!', "\"", '#', '$', '%', '&', "\'",
                     ')', '(', '*', '+', ',', '-', '.',
                     '/', '0', '1', '2', '3', '4', '5',
                     '6', '7', '8', '9', ':', ';', '<',
                     '=', '>', '?', '@', 'A', 'B', 'C',
                     'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K']

def getPhredToAccuDict():
    """
    Computes nucleotide sequencing accuracy of phred scores.
    0 to 1
    """
    return {k : 1 - pow(10, -i/10) for i, k in enumerate(phredSymbolList)}

def getPhredToErrorDict():
    """
    Computes nucleotide sequencing accuracy of phred scores.
    0 to 1
    """
    return {k : pow(10, -i/10) for i, k in enumerate(phredSymbolList)}

# Phred Level        A C T G
# 0-3    < 50 %      i n 1 ^
# 4-9    < 90 %      e h 2 &
# 10-19  < 99 %      o s 3 @
# 20-29  < 99.9 %    a c t g
# 30-42  < 99.99 %   A C T G

phredProbaseList = [
    4, 4, 4, 4, 3, 3, 3,
    3, 3, 3, 2, 2, 2, 2,
    2, 2, 2, 2, 2, 2, 1,
    1, 1, 1, 1, 1, 1, 1,
    1, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0]

baseToProbaseDict = {
    'A': ['A', 'a', 'o', 'e', 'i'],
    'T': ['T', 't', '3', '2', '1'],
    'C': ['C', 'c', 's', 'h', 'n'],
    'G': ['G', 'g', '@', '&', '^']
}
phredToProbaseIndicesDict = {i: k for i, k in zip(phredSymbolList, phredProbaseList)} 

def getProbaseFromBaseAndPhred(baseSeq, phredSeq):
    if len(baseSeq) != len(phredSeq):
        print("WARNING: Sequences of bases and their corresponding Phred scores differ in length.")
        return None
    
    probaseSeq = []
    for base, phredSymbol in zip(baseSeq, phredSeq):
        probaseIndex = phredToProbaseIndicesDict[phredSymbol]
        if base in ['A', 'T', 'C', 'G']:
            probaseSeq.append(baseToProbaseDict[base][probaseIndex])
        elif base == 'N':
            probaseSeq.append(base)
        else:
            "WARNING: Unknown base not in ATCGN."
    return ''.join(probaseSeq)

def getAvgPhredScore(phredSymbol1, phredSymbol2):
    return getNearestPhredFromAccu( (phredToAccuDict[phredSymbol1] + phredToAccuDict[phredSymbol2]) / 2 )

def getNearestPhredFromAccu(accuracyProb):
    phredQNumber = min( 42, round(-10 * log10(1 - accuracyProb)) )
    return phredSymbolList[phredQNumber]

def phredSeqToNpArray(phredSeq):
    return np.fromiter((phredToAccuDict[sym] for sym in phredSeq), float)

phredToAccuDict = getPhredToAccuDict();
phredToErrorDict = getPhredToErrorDict();