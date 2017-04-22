from manifestExtraction import genAmpliconDict
from difflib import SequenceMatcher
from csv import *

class AmpliconMatcherProbabilistic:
    def __init__(self, reference):
        self.ampDict1 = {}
        self.ampDict2 = {}
        self.refSeqs = []
        self.refSeqsTrunc = []
        self.ampDict1, self.ampDict2 = genAmpliconDict(reference)

        with open(reference) as csvFile:
            next(csvFile)
            self.refSeqs = [x[2] for x in list(reader(csvFile))]
            self.refSeqsTrunc = [x[:30] for x in self.refSeqs]
        print ("Done with init")

    def findAmpliconOld(self, read):
        """
        Deprecated
        Accepts a read of bases and returns the ID number of the amplicon associated with the read
        """
        ampID = str(self.ampDict1.get(read[:17], "000"))
        if ampID == "000":
            return str(self.ampDict2.get(read[7:24], "000"))
        else:
            return ampID

    def findAmplicon(self, read):
        """
        Accepts a read of bases and returns the ID number of the amplicon associated with the read
        """
        ampID = str(self.ampDict1.get(read[:17], "000"))
        if ampID == "000":
            ampID = str(self.ampDict2.get(read[7:24], "000"))

        # If the amplicon chosen are those which are easily mistaken, used cheem method to determine correct amplicon
        if ampID == "041" or ampID == "570":
            return self.findCorrect(read, self.refSeqs, [40, 569])
        elif ampID == "539" or ampID == "569":
            return self.findCorrect(read, self.refSeqs, [538, 568])
        elif ampID == "368" or ampID == "571":
            return self.findCorrect(read, self.refSeqs, [367, 570])
        elif ampID == "188" or ampID == "197":
            return self.findCorrect(read, self.refSeqs, [187, 196])
        elif ampID == "137" or ampID == "453":
            return self.findCorrect(read, self.refSeqs, [136, 452])
            """
            elif ampID == "000":

                possiblyCorrect = self.findCorrect(read[:30], self.refSeqsTrunc, range(571))
                if self.probDist(read, self.refSeqs[int(possiblyCorrect) - 1]) > 10000:
                    return possiblyCorrect
                else:
                    return "000"
            """
        else:
            return ampID

    def findCorrect(self, read, refSeqs, refsToCompare):
        return str(max([(self.probDist(read, self.refSeqs[x]), x) for x in refsToCompare])[1] + 1).rjust(3, "0")

    def probDist(self, read1, read2):
        return sum([1.2 ** x[2] for x in SequenceMatcher(None, read1, read2, autojunk = False).get_matching_blocks()])
