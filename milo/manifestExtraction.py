"""
Algorithms to extract selected data from manifest, and process them into useable forms.
"""

from csv import *
from itertools import *

def genTiledReads(inDir, outDir):
    """
    Combines overlapping amplicon tiles in the manifest into the contiguous sequenced region.
    Writes to 'IDs and Combined Alleles.csv
    """
    with open(inDir + "Manifest.csv") as csvFile:
        headerList = next(csvFile)
        csvList = list(reader(csvFile))
        orientedList = [[x[2], reverseComplement(x[3])] if x[1] == "+" else [x[3], reverseComplement(x[2])] for x in csvList]
        combinedList = (y[0][:y[0].rfind(y[1][:15])] + y[1] for y in orientedList)
        tileCounter = [(a, len(tuple(b))) for a, b in groupby([z[0][:-1] for z in csvList])]
        alleleList = []
        for c in tileCounter:
            for d in range(c[1]):
                if d == 0:
                    tempAllele = next(combinedList)
                else:
                    seqToAdd = next(combinedList)
                    tempAllele = tempAllele[:tempAllele.rfind(seqToAdd[:20])] + seqToAdd
            alleleList.append(tempAllele)
        outputList = list(zip([x[0] for x in tileCounter], alleleList))
        writeToCSV("IDs and Combined Alleles", outDir, outputList)

def genSignedReads(inDir, outDir):
    pass