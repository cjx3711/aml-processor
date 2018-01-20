"""
Tries to figure out which amplicon a sequence belongs to.
"""

class AmpliconMatcherHashSweep:

    def __init__(self, reference, kgramLength = 15, spacing = 15):
        """
        reference can be of two types, list or str.
        If list, it will use it directly, if str it will read the reference file
        """
        self.kgramLength = kgramLength
        self.spacing = spacing
        self.ampliconRefs = {}
        self.referenceCount = 0;
        
        # Config
        self.seqEndsWeight = 2
        self.seqMidWeight = 1

        self.noCount = 0
        self.badCount = 0
        self.ummCount = 0
        self.total = 0

        if type(reference) is list:
            self.generateReferenceFromList(reference)
        elif type(reference) is str:
            self.generateReferenceFromFile(reference)

        # Format of matchScore:
        # [
        #   [COUNT, ampID],
        #   ... 571 (number of references) times
        # ]
        self.resetMatchScore()

    def getWeight(self, cursorPos, readLen):
        """
        We want to weight the beginning of the sequence higher
        This is because the probe sequence is likely to be correct
        """
        if cursorPos > 30 or cursorPos < readLen - 30:
            return self.seqMidWeight
        else:
            return self.seqEndsWeight

    def calculateMatchScore(self, read):
        """
        Calculates the score for the read and returns the highest two
        possible reference amplicons
        """
        matches = []
        self.resetMatchScore()
        readCursor = 0 # The read cursor scans through the whole sequence one by one
        readLen = len(read)
        while readCursor < readLen - self.kgramLength + 1:
            kgram = read[readCursor:readCursor + self.kgramLength]
            if kgram in self.ampliconRefs:
                currentMatches = self.ampliconRefs[kgram]
                for match in currentMatches:
                    if (match in matches) == False: # This prevents multi counting.
                        ampID = match[0]
                        self.matchScore[ampID - 1][0] += self.getWeight(readCursor, readLen)
                        matches.append(match)

            readCursor += 1
        # Sort the scores so that the highest score is at the top
        self.matchScore.sort(reverse = True)
        return self.matchScore[0], self.matchScore[1]
        
    def findAmplicon(self, read):
        self.total += 1
        # largest1 = Amplicon with largest match score
        # largest2 = Amplicon with second largest match score
        # Format of largest1 & largest2: [COUNT, ampID]
        
        largest1, largest2 = self.calculateMatchScore(read)

        if ( largest1[0] < self.seqEndsWeight and largest2[0] == 0):
            self.noCount += 1
            return '000', None, 'nah'

        ampID = str(largest1[1]).rjust(3,'0')
        ampID2 = None
        matchType = 'kay'
        
        secondPercent = largest2[0] / largest1[0]
            
        if ( secondPercent > 0.5 ):
            self.badCount += 1
            ampID2 = str(largest2[1]).rjust(3,'0')
            matchType = 'bad'
            # umPerc = self.ummCount / self.total
            # badPerc = self.badCount / self.total
            # print("Umm, bad {0},{1} ({2},{3}) / {4}".format(self.ummCount, self.badCount, umPerc, badPerc, self.total))
        elif ( secondPercent > 0.3 ):
            self.ummCount += 1
            ampID2 = str(largest2[1]).rjust(3,'0')
            matchType = 'umm'
            # umPerc = self.ummCount / self.total
            # badPerc = self.badCount / self.total
            # print("Umm, bad {0},{1} ({2},{3}) / {4}".format(self.ummCount, self.badCount, umPerc, badPerc, self.total))
            
        return ampID, ampID2, matchType

    def generateReferenceFromFile(self, filename):
        with open(filename) as references:
            # Read through each line of the reference file and create the kgrams
            ampID = -1 # Skip the title line # 1-indexed
            for line in references:
                ampID += 1
                if ampID > 0:
                    self.referenceCount += 1
                    sequence = line.split(',')[3]
                    self.processSingleReferenceLine(sequence, ampID)

    def generateReferenceFromList(self, references):
        self.referenceCount = len(references)
        ampID = 0 # 1-indexed
        for sequence in references:
            ampID += 1
            self.processSingleReferenceLine(sequence, ampID)

    def processSingleReferenceLine(self, sequence, ampID):
        # Generate a bunch of kgrams as keys.
        for i in range(0,len(sequence) - self.kgramLength + 1, self.spacing):
            kgram = sequence[i:i+self.kgramLength]
            # Add each of the kgrams to the hash
            # A kgram can exist in multiple amplicons, so we store the amplicons it belongs to
            if (kgram in self.ampliconRefs) == False:
                self.ampliconRefs[kgram] = []
            self.ampliconRefs[kgram].append((ampID,i))

    def getReferenceCount(self):
        return self.referenceCount
        
    def resetMatchScore(self):
        self.matchScore = [[0,x + 1] for x in range(self.referenceCount)]
