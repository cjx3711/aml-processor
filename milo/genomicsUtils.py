"""
genomicsUtils | Methods for common routines in genome analysis
"""
from math import pow
from csv import *
import os

complement = {"A" : "T", "C" : "G", "G" : "C", "T" : "A", " " : " "}
phredList = [
                "K", "J", "I", "H", "G", "F", "E", 
                "D", "C", "B", "A", "@", "?", ">", 
                "=", "<", ";", ":", "9", "8", "7", 
                "6", "5","4", "3", "2", "1", "0", 
                "/", ".", "-", ",", "+", "*", "(", 
                ")", "\'", "&", "%", "$", "#", "\"", "!"
                ]

def reverseComplement(bases):
    """
    Reverse complements an iterable of nucleotides.
    Non ATCG nucleotides are converted to N
    """
    return ''.join(reversed([complement.get(x, "N") for x in bases]))

def writeToCSV(fileName, outdir, targetList):
    """
    Creates filename.csv in outdir, or overrides existing file if any.
    Writes *targetList as rows
    """
    with open(fileName + ".csv", "w+", newline = "") as outFile:
        outWriter = writer(outFile)
        outWriter.writerows(targetList)
        outFile.close()

def grouper(iterable, n):
    """
    Generator which returns a zip object containing the next n values in iterable.
    """
    args = [iter(iterable)] * n
    return zip(*args)

def nthAndKthLetter(targetString, targetLetter, n, k):
    """
    Finds the nth and kth occurence of targetLetter in targetString.
    """
    lastIndex, nIndex = 0, 0
    for i in range(k):
        if i == n: nIndex = lastIndex
        lastIndex = targetString.index(targetLetter, lastIndex) + 1
    return (nIndex, lastIndex - 1)

def getPhredQualityDict():
    """
    Computes nucleotide sequencing accuracy of phred scores.
    0 to 1
    """
    return {x : 1 - pow(10, -y/10) for x, y in zip(phredList, range(42, -1, -1))}

def getPhredErrorDict():
    """
    Computes nucleotide sequencing accuracy of phred scores.
    0 to 1
    """
    return {x : pow(10, -y/10) for x, y in zip(phredList, range(42, -1, -1))}
    
def simpleDistance(a, b):
    """
    Deprecated: Replaced by fastcomp
    Checks the distance between two strings.
    Outputs 0, 1 or 2. If it's 2, it means it's more than 2
    """
    if ( a == b ):
        return 0
    if (len(a) != len(b)):
        return 2
    cursorA = 0
    cursorB = 0
    differences = 0
    
    while ( cursorA < len(a) and cursorB < len(b) ):
        if ( a[cursorA] != b[cursorB] ):
            differences += 1
            if ( differences >= 2 ):
                break;
        cursorA += 1
        cursorB += 1
    
    return differences

def extractChromosomeNumber(ampliconName):
    nameParts = ampliconName.split('.')
    for namePart in nameParts:
        if namePart.startswith('chr'):
            return namePart[3:]
    return ''