"""
Replaces ReadPairer.
Makes use of a more probabilistic method to determine the read pairing.
JJ please expound on this.
"""

from collections import deque
from phredUtils import *
from DataTypes import *

class PhredBasedReadPairer:
    def __init__(self, configFile = 'config.json'):
        self.MATCH_SCORE = 1
        self.MISMATCH_PENALTY = 5        
        self.SCORE_THRESHOLD = 10

        # WE SHOULD MOVE THIS OUT. Check the first line of FASTQ1 and FASTQ2, and assert if they are the same.
        self.READ_LENGTH = 151
        self.SUPERPOSITION_INDEX = self.READ_LENGTH - self.SCORE_THRESHOLD
        

    def pairRead(self, leftSeq, leftPhred, rightSeq, rightPhred):
        """
        Takes two sequences (bases and respective Phred scores) which have already been reverse complemented
        and merges them into a paired read. Returns a PairedRead namedtuple.
        """
        # Get the smallest overlapping substring of bases which can fulfil the score threshold, and its Phred qualities
        leftSubstr = deque(leftSeq[-self.SCORE_THRESHOLD:])
        leftSubstrQual = deque(leftPhred[-self.SCORE_THRESHOLD:])
        rightSubstr = deque(rightSeq[:self.SCORE_THRESHOLD])
        rightSubstrQual = deque(rightPhred[:self.SCORE_THRESHOLD])

        dequeLen = self.SCORE_THRESHOLD

        # Iteratively grow the overlapping region by sliding the left and right reads over each other
        for i in range(0, self.SUPERPOSITION_INDEX):
            score = 0
            for j in range(dequeLen):
                score += self.calcScore(leftSubstr[j], leftSubstrQual[j], rightSubstr[j], rightSubstrQual[j])
            # If score meets threshold, merge the overlapping region, then append the non-overlapping ends and return
            if score >= self.SCORE_THRESHOLD:
                mergedRead = self.mergeOverlap(leftSubstr, leftSubstrQual, rightSubstr, rightSubstrQual)
                return PairedRead(leftSeq[:self.SUPERPOSITION_INDEX - i] + mergedRead[0] + rightSeq[self.SCORE_THRESHOLD + i:],
                        leftPhred[:self.SUPERPOSITION_INDEX - i] + mergedRead[1] + rightPhred[self.SCORE_THRESHOLD + i:],
                        True)
            else:
                leftIndexToAdd = -self.SCORE_THRESHOLD - i - 1
                rightIndexToAdd = self.SCORE_THRESHOLD + i

                leftSubstr.appendleft(leftSeq[leftIndexToAdd])
                leftSubstrQual.appendleft(leftPhred[leftIndexToAdd])
                rightSubstr.append(rightSeq[rightIndexToAdd])
                rightSubstrQual.append(rightPhred[rightIndexToAdd])
                dequeLen += 1

        # When the sliding crosses the halfway (superposition) mark, we have to shrink the overlapping region instead
        for _ in range(0, self.SUPERPOSITION_INDEX):
            for j in range(dequeLen):
                score += self.calcScore(leftSubstr[j], leftSubstrQual[j], rightSubstr[j], rightSubstrQual[j])
            if score >= self.SCORE_THRESHOLD:
                # Note that there are no overlapping ends to append in this case
                baseSeq, phredSeq = self.mergeOverlap(leftSubstr, leftSubstrQual, rightSubstr, rightSubstrQual)
                return PairedRead(baseSeq, phredSeq, True) 
            else:
                leftSubstr.pop()
                leftSubstrQual.pop()
                rightSubstr.popleft()
                rightSubstrQual.popleft()
                dequeLen -= 1
            
        return PairedRead(leftSeq + '|' + rightSeq, leftPhred + '|' + rightPhred, False)

    def calcScore(self, base1, phred1, base2, phred2):
        """
        Gives a penalty/reward based on mis/match of two bases, weighted by quality score
        """
        scoreWeight = min(phredToAccuDict[phred1], phredToAccuDict[phred2]) ** 2
        if base1 == base2:
            return self.MATCH_SCORE * scoreWeight
        elif base1 == 'N' or base2 == 'N':
            return 0
        else:
            return - self.MISMATCH_PENALTY * scoreWeight

    def mergeOverlap(self, leftSubstr, leftSubstrPhred, rightSubstr, rightSubstrQual):
        mergedBases = []
        mergedQual = []
        for i in range(len(leftSubstr)):
            # If both bases are the same, then just pick the base and average the qualities
            if leftSubstr[i] == rightSubstr[i]:
                mergedBases.append(leftSubstr[i])
                mergedQual.append(getAvgPhredScore(leftSubstrPhred[i], rightSubstrQual[i]))
            # If bases are different, but left base is more accurate, pick it as the base but dilute the quality score
            elif phredToErrorDict[leftSubstrPhred[i]] < phredToErrorDict[rightSubstrQual[i]]:
                # An approximation of the better base quality, ignoring inferior base's weight on error distribution
                probBetterBaseIsCorrect = (phredToErrorDict[rightSubstrQual[i]] / 3 + phredToAccuDict[leftSubstrPhred[i]]) / 2
                mergedQual.append(getNearestPhredFromAccu(probBetterBaseIsCorrect))
                mergedBases.append(leftSubstr[i])
            # If the right base is more accurate...
            else:
                probBetterBaseIsCorrect = (phredToErrorDict[leftSubstrPhred[i]] / 3 + phredToAccuDict[rightSubstrQual[i]]) / 2
                mergedQual.append(getNearestPhredFromAccu(probBetterBaseIsCorrect))
                mergedBases.append(rightSubstr[i])
        return ''.join(mergedBases), ''.join(mergedQual)