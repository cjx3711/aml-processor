"""
genomicsUtils | Methods for common routines in genome analysis
"""
from math import pow

complement = {"A" : "T", "C" : "G", "G" : "C", "T" : "A"}
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
    """
    return {x : 1 - pow(10, -y/10) for x, y in zip(phredList, range(42, -1, -1))}

# Temp unit for testing the testing framework
def IsOdd(n):
    return n % 2 == 1