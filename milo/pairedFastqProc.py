from Bio.pairwise2 import *
from j3xUtils import *
from concurrent.futures import *
from itertools import *

numThreads = 12
baseScore = {"A": 1, "T": 1, "C": 1, "G": 1, "N": 0}
readLength = 151
scoreThreshold = 10
rangeEnd = readLength - scoreThreshold
rangeStart = -rangeEnd

def mergeUnpaired(left, right, lQuality, rQuality):
    #Align unpaired reads
    for x in range(rangeStart, rangeEnd):
        # x is the relative position of the left read to the right, as it slides rightward
        # lRange and rRange contain the (start, end) indices of the overlap in the left and right reads
        lRange = (max(0, -x), min(readLength, -x + readLength))
        rRange = (max(0, x), min(readLength, x + readLength))
        # Generates 2D tuples of paired (left base, right base), and (left quality, right quality)
        overlapPairs = tuple(zip(left[lRange[0]:lRange[1]], right[rRange[0]:rRange[1]]))
        overlapQuality = tuple(zip(lQuality[lRange[0]:lRange[1]], rQuality[rRange[0]:rRange[1]]))
        overlapLength = lRange[1] - lRange[0]
        # If score in overlapping region is high, combine left and right reads
        if calcScore(left, right, overlapPairs, overlapLength) > scoreThreshold:
            # mergedSeq = ((bases, quality), collisions)
            mergedSeq = mergeOverlap(overlapPairs, overlapQuality)
            # Append the non-overlapping bases to the left and right of the merged sequence
            allBases = chain(left[0:lRange[0]], mergedSeq[0][0], right[rRange[1]:readLength])
            allQuality = chain(lQuality[0:lRange[0]], mergedSeq[0][1], rQuality[rRange[1]:readLength])
            # Convert to j3x format
            j3xBases = "".join("_" if x == "N" else x for x in allBases)
            j3xQuality = "".join(simplifyQuality(qualityDict[x]) for x in allQuality)
            return j3xBases, j3xQuality, str(mergedSeq[1])
    # If alignment fails, return both reads separated with a space
    return " ".join((left, right)), " ".join((lQuality, rQuality)), "N/A"

def calcScore(left, right, overlapPairs, overlapLength):
    score, count = 0, 0
    for x, y in overlapPairs:
        # Rewards exact matches, heavily penalizes differences, ignores when quality is too low
        if x == "N" or y == "N":
            pass
        elif x == y:
            score += 1
        elif x != y:
            score -= 5
        count += 1

        # If score is too low to possibly pass required score, break
        if (score + overlapLength - count) <= scoreThreshold:
            # Note: Low quality scores are not correctly calculated because of this
            break
    return score

def mergeOverlap(overlapPairs, overlapQuality):
    # Collisions are the number of contradictions during a merge, defined as bases whose quality scores do not diverge significantly
    collisions = [0]
    def pickBetter(bases, quality):
        """
        Helper function to pick the base with better quality
        """
        if bases[0] != bases[1]:
            # If the bases are different
            if errorDict[quality[0]] > errorDict[quality[1]] * 5:
                # If the right read is at least five times less likely to have an error, return it and its quality
                return bases[1], quality[1]
            elif errorDict[quality[1]] > errorDict[quality[0]] * 5:
                # Otherwise, just return the other read
                pass
            else:
                # If one is not a significantly better choice than the other, increment collision counter.
                collisions[0] += 1
        return bases[0], quality[0]
    return tuple("".join(y) for y in zip(*(pickBetter(*x) for x in zip(overlapPairs, overlapQuality)))), collisions[0]

def alignAndMerge(left, right):
    # Merges (left sequence, reverseComplement(right sequence), left quality, reversed(right quality))
    basesQualityCollisions = mergeUnpaired(left[1][:-1], reverseComplement(right[1][:-1]), left[3][:-1], right[3][:-1][::-1])
    # Creates new ID sequence by retrieving only the coordinates from the existing FASTQ read ID, and then appending the number of collisions
    coordIndices = nthAndKthLetter(left[0], ":", 5, 7)
    sequenceID = left[0][coordIndices[0]: coordIndices[1] - 2]
    newID = ["".join(("C:", basesQualityCollisions[2], ", ", sequenceID))]

    # Adds the sequence and quality, and returns
    newID.extend(basesQualityCollisions[:2])
    return "\n".join(newID)

def pairToJ3X(fq1, fq2, inDir, outDir):
        with open(inDir + fq1) as fq1File, open(inDir + fq2) as fq2File:
            with open(outDir + fq1[:(fq1.find("_R1_"))] + "PAIRED.j3x", "w+", newline = "") as outFile:
                # Creates iterators which deliver the 4 lines of each FASTQ read as a zip (ID, Sequence, Blank, Quality)
                fq1Iter, fq2Iter = grouper(fq1File, 4), grouper(fq2File, 4)
                with ProcessPoolExecutor(numThreads) as processManager:
                    # Calls alignAndMerge(FASTQ1's (ID, Sequence, Blank, Quality), FASTQ2's (ID, Sequence, Blank, Quality))
                    for x in processManager.map(alignAndMerge, fq1Iter, fq2Iter, chunksize = 250):
                        outFile.write(x)
                        outFile.write("\n\n")
                    outFile.close()