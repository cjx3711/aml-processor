class AmpliconMatcherKgram:

    def __init__(self, manifest, kgramLength = 15, spacing = 15):
        self.kgramLength = kgramLength
        self.spacing = spacing
        self.indelHash = {}
        self.referenceCount = 0;

        self.noCount = 0
        self.badCount = 0
        self.ummCount = 0
        self.total = 0
        self.generateManifest(manifest)


    def generateManifest(self, filename):
        with open(filename, "r", newline = "") as references:

            # Read through each line of the reference file and create the kgrams
            lineno = -1
            for line in references:
                if lineno >= 0:
                    self.referenceCount += 1
                    sequence = line.split(',')[2]
                    for i in range(0,len(sequence) - self.kgramLength, self.spacing):
                        kgram = sequence[i:i+self.kgramLength]
                        # Add each of the kgrams to the hash
                        if (kgram in self.indelHash) == False:
                            self.indelHash[kgram] = []
                        self.indelHash[kgram].append((lineno,i))
                lineno += 1
        self.matchCounts = [[0,x] for x in range(self.referenceCount)]

    def resetMatchCounts(self):
        i = 0
        for mc in self.matchCounts:
            mc[0] = 0
            mc[1] = i
            i += 1

    def findAmpliconKgram(self, read):
        self.total += 1
        matches = []
        self.resetMatchCounts()
        i = 0
        while i < len(read) - self.kgramLength + 1:
            kgram = read[i:i+self.kgramLength]
            if kgram in self.indelHash:
                currentMatches = self.indelHash[kgram]
                i += self.spacing - 1
                for match in currentMatches:
                    self.matchCounts[match[0]][0] += 1
                matches.extend(currentMatches)
            i += 1
        self.matchCounts.sort(reverse = True)
        # largest1, largest2 = getLargestTwo(matchCounts)
        secondPercent = 0 # Percentage of the second largest
        if ( self.matchCounts[0][0] == 0 and self.matchCounts[1][0] == 0):
            self.noCount += 1
            return '000'

        if ( self.matchCounts[0][0] != 0 ):
            secondPercent = self.matchCounts[1][0] / self.matchCounts[0][0]
        if ( secondPercent > 0.6 ):
            self.badCount += 1
        elif ( secondPercent > 0.3 ):
            self.ummCount += 1

        return str(self.matchCounts[0][1] + 1).rjust(3,'0')
