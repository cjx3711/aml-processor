from j3xUtils import *
from manifestExtraction import *
from concurrent.futures import *
from itertools import *
from difflib import *

numThreads = 12
baseScore = {"A": 1, "T": 1, "C": 1, "G": 1, "N": 0}
readLength = 151
scoreThreshold = 10
rangeEnd = readLength - scoreThreshold
rangeStart = -rangeEnd
ampDict1, ampDict2 = genAmpliconDict("references/")

def mergeUnpaired(left, right, lQuality, rQuality):
    """
    Note: ulnernable to early terminaton upon encountering random repeats or homopolymer sequences
    Merges two unpaired reads
    Quality is not taken into consideration for tha calculations, it's only used as a return.
    Returns mergedSequence, qualityScores, numberOfCollisions
    """
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
        if calcScore(overlapPairs, overlapLength) > scoreThreshold:
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
    return " ".join((left, right)), " ".join((lQuality, rQuality)), "?"

def calcScore(overlapPairs, overlapLength):
    """
    Checks if the overlap is legit.
    Returns a nice score for the overlap
    """
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

def findAmplicon(read):
    """
    Accepts a read of bases and returns the ID number of the amplicon associated with the read
    """
    ampID = str(ampDict1.get(read[:17], "000"))
    if ampID == "000":
        return str(ampDict2.get(read[7:24], "000"))
    else:
        return ampID

with open("references/Manifest.csv") as csvFile:
        next(csvFile)
        refSeqs = [x[2] for x in list(reader(csvFile))]

def findAmpliconLevenshtein(read):
    """
    Accepts a read of bases and returns the ID number of the amplicon associated with the read
    """
    ampID = str(ampDict1.get(read[:17], "000"))
    if ampID == "000":
        ampID = str(ampDict2.get(read[7:24], "000"))

    # If the amplicon chosen are those which are easily mistaken, used Levenshtein to determine correct amplicon
    # CAUTION: Vulnerable to large insertions
    if ampID == "041" or ampID == "570":
        return findCorrect(read, refSeqs, 41, 570)
    elif ampID == "539" or ampID == "569":
        return findCorrect(read, refSeqs, 539, 569)
    elif ampID == "368" or ampID == "571":
        return findCorrect(read, refSeqs, 368, 571)
    elif ampID == "188" or ampID == "197":
        return findCorrect(read, refSeqs, 188, 197)
    elif ampID == "137" or ampID == "453":
        return findCorrect(read, refSeqs, 137, 453)
    else:
        return ampID

def findCorrect(read, refSeqs, i, j):
    i -= 1
    j -= 1
    return str(min((distance(read, refSeqs[i]), i),(distance(read, refSeqs[j]), j))[1] + 1).rjust(3, "0")

def probDist(read1, read2):
    return sum([1.5 ** x[2] for x in SequenceMatcher(None, read1, read2, autojunk = False).get_matching_blocks()])

def alignAndMerge(left, right):
    # Merges (left sequence, reverseComplement(right sequence), left quality, reversed(right quality))
    basesQualityCollisions = mergeUnpaired(left[1][:-1], reverseComplement(right[1][:-1]), left[3][:-1], right[3][:-1][::-1])
    # Retrieves the coordinates from the existing FASTQ read ID
    coordIndices = nthAndKthLetter(left[0], ":", 5, 7)
    sequenceID = left[0][coordIndices[0]: coordIndices[1] - 2]
    # Checks which amplicon a read belongs to
    ampID = findAmpliconLevenshtein(basesQualityCollisions[0])
    # Joins amplicon number, collision number, and coordinate as new ID
    newID = ["".join(("ID:", ampID, ", C:", basesQualityCollisions[2], ", ", sequenceID))]
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