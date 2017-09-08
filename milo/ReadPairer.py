from j3xUtils import *
from manifestExtraction import *
from concurrent.futures import *
from DataTypes import *
from itertools import *
from AmpliconMatcherHashSweep import *
import json

class ReadPairer:
    def __init__(self, configFile = 'config.json', referenceFile = 'references/Manifest.csv'):
        self.alignByMaxima = False # Whether or not to robustly align by detecting global maxima for overlap scores
        self.matchScore = 1
        self.mismatchPenalty = 5
        with open(configFile) as config_file: 
            config_data = json.load(config_file)
            if ( 'readPairer_alignByMaxima' in config_data ):
                self.alignByMaxima = config_data['readPairer_alignByMaxima']
            if ( 'readPairer_matchScore' in config_data ):
                self.matchScore = abs(config_data['readPairer_matchScore'])
            if ( 'readPairer_mismatchPenalty' in config_data ):
                self.mismatchPenalty = abs(config_data['readPairer_mismatchPenalty'])
        
        self.readLength = 151
        self.scoreThreshold = 5
        self.rangeEnd = self.readLength - self.scoreThreshold
        self.rangeStart = -self.rangeEnd
        self.ampliconMatcher = AmpliconMatcherHashSweep(referenceFile)

    def getReferenceCount(self):
        return self.ampliconMatcher.getReferenceCount()

    def mergeUnpaired(self, left, right, lQuality, rQuality, alignByMaxima):
        newScore = 0 # Overlap score in current iteration
        maxScore = 0 # Local maxima for overlap similarity
        bestMatch = [] # Stores the overlap coordinates of local maxima of overlap score [lRange, rRange]
        # Slides left read rightward
        for x in range(self.rangeStart, self.rangeEnd):
            lRange = (max(0, -x), min(self.readLength, -x + self.readLength)) # The indices of the overlap zone on the left read
            rRange = (max(0, x), min(self.readLength, x + self.readLength))
            newScore = self.calcScore(self.getOverlapPairs(left, right, lRange, rRange), lRange[1] - lRange[0], maxScore, alignByMaxima)
            
            if alignByMaxima:
                if newScore > maxScore: # If we encounter a new local maxima, record the score and its overlap coordinates
                    maxScore = newScore
                    bestMatch = [lRange, rRange]
            elif newScore > self.scoreThreshold: # If we aren't checking for the maxima, then whenever we exceed the score threshold, merge and return
                return self.alignCoordsToJ3x(left, right, lQuality, rQuality, lRange, rRange)
        
        if alignByMaxima and maxScore > 0: # If we are aligning by maxima, and there exists a global maxima better than not pairing at all
            lRange, rRange = bestMatch
            return self.alignCoordsToJ3x(left, right, lQuality, rQuality, lRange, rRange)
        else: # If we cannot pair properly, return both reads separated with a space
            return " ".join((left, right)), " ".join((lQuality, rQuality)), "?"

    def calcScore(self, overlapPairs, overlapLength, maxScore, alignByMaxima):
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
                score += self.matchScore
            elif x != y:
                score -= self.mismatchPenalty
            count += 1

            # If the remaining bases can never give a score higher than the maxima/threshold, break
            if alignByMaxima:
                if (score + overlapLength - count) <= maxScore:
                    break
            elif (score + overlapLength - count) <= self.scoreThreshold:
                break
        return score

    def getOverlapPairs(self, left, right, lRange, rRange):
        # Takes the overlapping regions from left and right, and returns tuple of opposing pairs
        return zip(left[lRange[0]:lRange[1]], right[rRange[0]:rRange[1]])

    def alignCoordsToJ3x(self, left, right, lQuality, rQuality, lRange, rRange):
        # Get the merged reads and quality sequences
        mergedSeq, collisions = self.mergeOverlap(self.getOverlapPairs(left, right, lRange, rRange), self.getOverlapPairs(lQuality, rQuality, lRange, rRange))
        # Paste the non-overlapping parts at the front and back
        allBases = chain(left[0:lRange[0]], mergedSeq[0], right[rRange[1]:self.readLength])
        allQuality = chain(lQuality[0:lRange[0]], mergedSeq[1], rQuality[rRange[1]:self.readLength])
        # Convert to j3x format
        j3xBases = "".join("_" if x == "N" else x for x in allBases)
        j3xQuality = "".join(simplifyQuality(qualityDict[x]) for x in allQuality)
        return j3xBases, j3xQuality, str(collisions)

    def mergeOverlap(self, overlapPairs, overlapQuality):
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

    def alignAndMerge(self, left, right):
        bases, quality, collisions = self.mergeUnpaired(left[1].rstrip(), reverseComplement(right[1].rstrip()), left[3].rstrip(), right[3].rstrip()[::-1], self.alignByMaxima)
        
        failedToPair = 1 if collisions == '?' else 0
        
        # # Retrieves the coordinates from the existing FASTQ read ID
        # coordIndices = nthAndKthLetter(left[0], ":", 5, 7)
        # sequenceID = left[0][coordIndices[0]: coordIndices[1] - 2]
        # Checks which amplicon a read belongs to
        ampID, ampID2, matchType = self.ampliconMatcher.findAmplicon(bases)
        
        
        # Joins amplicon number, collision number, and coordinate as new ID, and appends bases and quality
        IDPart = ''
        if ampID2 != None:
            IDPart = 'TL:{0}/{1}'.format(ampID, ampID2)
        else:
            IDPart = 'ID:{0}'.format(ampID)
            
        readData = ", ".join((IDPart, 'C:'+collisions))
        return AlignedAndMerged(failedToPair, matchType, readData, bases, quality)

        # otherStats = (failedToPair, matchType)
        # return (otherStats, ["".join(("ID:", ampID, ", C:", collisions, ", ", sequenceID)), bases, quality])
            
