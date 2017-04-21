"""
Algorithms to extract selected data from manifest, and process them into useable forms.
"""
from itertools import *
from genomicsUtils import *

def genTiledReads(inDir, outDir):
    """
    Combines overlapping amplicon tiles in the manifest into the contiguous sequenced region. Writes result to CSV.
    """
    with open(inDir + "Manifest.csv") as csvFile:
        # Puts column labels in headerList and rest of CSV file in csvList
        headerList = next(csvFile)
        csvList = list(reader(csvFile))

        # Reverse complements either the left or right read depending on strand orientation, puts corrected [leftRead, rightRead] in orientedList
        orientedList = [[x[2], reverseComplement(x[3])] if x[1] == "+" else [x[3], reverseComplement(x[2])] for x in csvList]
        # Merges both reads by adding leftRead[:rightRead starts] + rightRead, puts in combinedList
        combinedList = [left[:left.rfind(right[:15])] + right for left, right in orientedList]
        # Groups amplicons of same target, puts [target, number of amplicons covering target] into tileCounter
        tileCounter = [(target, len(tuple(group))) for target, group in groupby([amp[0][:amp[0].rfind("_")] for amp in csvList])]
        targetList = []
        for target, numAmps in tileCounter:
            for i in range(numAmps):
                # Repeat this for the total number of tiles in each target
                if i == 0:
                    # When processing a new target, reset our temporary target sequence to tile_1 of the target
                    tempTargetSeq = next(combinedList)
                else:
                    # Else, merge the next tile with our temporary target sequence
                    seqToAdd = next(combinedList)
                    tempTargetSeq = tempTargetSeq[:tempTargetSeq.rfind(seqToAdd[:20])] + seqToAdd
            # When all tiles in target have been merged, add to targetList
            targetList.append(tempAllele)

        # Writes [target, merged target sequence] to CSV
        outputList = list(zip([x[0] for x in tileCounter], targetList))
        writeToCSV("IDs and Combined Alleles", outDir, outputList)

def genSignedReads(inDir, outDir):
    """
    Merges left and right reference reads into signed amplicon reference. Writes result to CSV.
    """
    with open(inDir + "Manifest.csv") as csvFile:
        # Puts column labels in headerList and rest of CSV file in csvList
        headerList = next(csvFile)
        csvList = list(reader(csvFile))

        # Puts [leftRead, reverseComplement(rightRead)] in diffStrands
        diffStrands = [[x[2], reverseComplement(x[3])] for x in csvList]
        # Merges both reads by adding leftRead[:rightRead starts] + rightRead
        combinedList = [left[:left.rfind(right[:15])] + right for left, right in diffStrands]

        # Writes [leftRead, rightRead, mergedAmplicon] to CSV
        outputList = list(zip(diffStrands, combinedList))
        writeToCSV("Paired Signed Reads", outDir, outputList)

def genAmpliconDict(inDir):
    """
    CAUTION
    Generates a dictionary of amplicon hashes. For testing purposes only.
    """
    with open(inDir + "Manifest.csv") as csvFile: